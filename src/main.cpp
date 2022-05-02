

#include<iostream>
#include<algorithm>
#include<string>
#include<sstream>
#include<cmath>
#include<vector>
#include<fstream>
#include<cassert>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


using namespace std;

double Mean(vector<double> v){
    double sum = 0;
    for(int n=0;n<v.size();n+=1){
        sum += v[n];
    }
    double mean = sum / v.size();
    return mean;
}

vector<int> GetMatchPos(vector<vector<int>> fingerings){
	assert(fingerings.size()>1);
	for(int k=1;k<fingerings.size();k+=1){
		assert(fingerings[k].size()==fingerings[0].size());
	}//endfor k

	vector<int> matchPos;
	for(int n=0;n<fingerings[0].size();n+=1){
		bool allMatch=true;
		for(int k=1;k<fingerings.size();k+=1){
			if(fingerings[k][n]!=fingerings[0][n]){
				allMatch=false;
				break;
			}//endif
		}//endfor k
		if(allMatch){matchPos.push_back(n);}
	}//endfor n

	return matchPos;
}//end GetMatchPos


#include <typeinfo>

double MultiGTError(vector<vector<int>> finsGT,vector<int> finEst,double substCost=1,double softSwitchCost=0.1,double hardSwitchCost=1000){
	assert(finsGT.size()>0);
	for(int k=1;k<finsGT.size();k+=1){
		assert(finsGT[k].size()==finsGT[0].size());
	}//endfor k
	assert(finEst.size()==finsGT[0].size());


	int nGT=finsGT.size();

	vector<vector<int> > seqEst(2);//[0]=RH,[1]=LH
	vector<vector<vector<int> > > seqGT(2);//[0]=RH,[1]=LH
	for(int h=0;h<2;h+=1){seqGT[h].resize(nGT);}

	for(int n=0;n<finEst.size();n+=1){
		int fn=finEst[n];
		int hand=0;
		if(fn<0){hand=1;}//endif LH
		seqEst[hand].push_back(fn);
		for(int k=0;k<nGT;k+=1){
			fn=finsGT[k][n];
			seqGT[hand][k].push_back(fn);
		}//endfor k
	}//endfor n


	vector<double> cumuCost(2);

	for(int h=0;h<2;h+=1){
		int len=seqEst[h].size();
		if(len != 0){
            vector<double> cost(nGT);
            vector<vector<int> > amin(len);
            vector<int> optPath(len);//association to GT;

            for(int n=0;n<len;n+=1){
                amin[n].resize(nGT);
                if(n==0){
                    cost.assign(nGT,0);
                    for(int z=0;z<nGT;z+=1){
                        if(seqEst[h][n]!=seqGT[h][z][n]){cost[z]+=substCost;}
                    }//endfor z
                    continue;
                }//endif

                vector<double> preCost(cost);
                double tmpCost;
                for(int z=0;z<nGT;z+=1){
                    cost[z]=preCost[z];
                    amin[n][z]=z;
                    for(int zp=0;zp<nGT;zp+=1){
                        if(zp==z){continue;}
                        tmpCost=preCost[zp]+((seqGT[h][zp][n-1]==seqGT[h][z][n-1])? softSwitchCost:hardSwitchCost);
                        if(tmpCost<cost[z]){
                            cost[z]=tmpCost;
                            amin[n][z]=zp;
                        }//endif
                    }//endfor zp
                    if(seqEst[h][n]!=seqGT[h][z][n]){cost[z]+=substCost;}
                }//endfor z
            }//endfor n
            optPath[len-1]=0;
            for(int z=0;z<nGT;z+=1){
                if(cost[z]<cost[optPath[len-1]]){optPath[len-1]=z;}
            }//endfor z
            for(int n=len-2;n>=0;n-=1){
                optPath[n]=amin[n+1][optPath[n+1]];
            }//endfor n
            cumuCost[h]=cost[optPath[len-1]];
        }
        else{
                cumuCost[h] = 0.0;
        }
    }//endfor h
	return cumuCost[0]+cumuCost[1];
}//end MultiGTError


double AveragePairwiseMatchRate(vector<vector<int>> finsGT,vector<int> finEst){
	assert(finsGT.size()>0);
	for(int k=1;k<finsGT.size();k+=1){
		assert(finsGT[k].size()==finsGT[0].size());
	}//endfor k
	assert(finEst.size()==finsGT[0].size());

	int nGT=finsGT.size();
	int nNotes=finEst.size();

	vector<double> matchRates(nGT);

	for(int k=0;k<finsGT.size();k+=1){
		vector<vector<int>> pair;
		pair.push_back(finsGT[k]);
		pair.push_back(finEst);
		vector<int> matchPos=GetMatchPos(pair);
		matchRates[k]=double(matchPos.size())/nNotes;
	}//endfor k

	return Mean(matchRates);
}//end AveragePairwiseMatchRate

double hmr(vector<vector<int>> finsGT,vector<int> finEst){
    return (finEst.size()-MultiGTError(finsGT,finEst,1,10000,10000))/double(finEst.size());
}

double smr(vector<vector<int>> finsGT,vector<int> finEst){
    return (finEst.size()-MultiGTError(finsGT,finEst,1,0,0))/double(finEst.size());
}

double rmr(vector<vector<int>> finsGT,vector<int> finEst){
    return (finEst.size()-MultiGTError(finsGT,finEst,1,1,10000))/double(finEst.size());
}



namespace py = pybind11;

PYBIND11_MODULE(multi_fingering_eval, m) {
    // optional module docstring
    m.doc() = "nakamura metrics pybind";

    // define add function
    m.def("gmr", &AveragePairwiseMatchRate, "general match rate");
    m.def("hmr", &hmr, "highest match rate");
    m.def("smr", &smr, "soft match rate");
    m.def("rmr", &rmr, "recomination match rate");
}

