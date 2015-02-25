/*
 * =====================================================================================
 *
 *       Filename: auxi.cpp
 *    Description: 
 *        Version: 1.0
 *        Created: 08/22/2014 03:57:23 PM
 *  Last Modified: 01/29/2015 03:39:42 PM
 *       Compiler: gcc
 *
 *         Author: Wei Zheng
 *          Email: intuitionofmind@gmail.com
 *   Organization: Institute for Advanced Study, Tsinghua University
 *
 * =====================================================================================
 */

#include"head.h"
#include"para.h"

int num_site(string s){
    if("PBC"==s){
        return num_dx*num_dy;
    }
    else{
        if("CBC"==s){
            return (num_dx*num_dy+num_dy);
        }
        else{
            return -1;
        }
    }
}

double unit_prob(){
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<double> dis(0, 1);
    return dis(generator);
}

int search(string s, int x, int y){  //s is the indicator to the boundary condition: periodic or cylinderical
    int index;
    if("PBC"==s){
        index=y*num_dx+x;
        return index;
    }
    else{
        if("CBC"==s){
            index=y*(num_dx+1)+x;
            return index;
        }
        else{
            return -1;
        }
    }
}

int initialize_spin(int len, int* Spin){
    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<int> dis(0, 1);
    for(int i=0; i<len; i++){
        Spin[i]=(dis(generator))*2-1;
    }
    return 1;
}

bool locate_bond(string s, int b, int* Coor){
    if(b<(num_dx*num_dy)){
        Coor[0]=b%num_dx;
        Coor[1]=b/num_dx;
        if("PBC"==s){
            Coor[2]=(Coor[0]+1)%num_dx;
        }
        else{
            if("CBC"==s){
                Coor[2]=Coor[0]+1;
            }
            else{
                return false;
            }
        }
        Coor[3]=Coor[1];
        return true;
    }
    else{
        int bb=b-(num_dx*num_dy);
        if("PBC"==s){
            Coor[0]=bb%num_dx;
            Coor[1]=bb/num_dx;
            Coor[2]=Coor[0];
            Coor[3]=(Coor[1]+1)%num_dy;
            return true;
        }
        else{
            if("CBC"==s){
                Coor[0]=bb%(num_dx+1);
                Coor[1]=bb/(num_dx+1);
                Coor[2]=Coor[0];
                Coor[3]=(Coor[1]+1)%num_dy;
                return true;
            }
            else{
                return false;
            }
        }
    }
}

double check_bond(string bc, int b){
    int* Coor=new int[num_leg];
    locate_bond(bc, b, Coor);
    if((b<num_dx*num_dy)&&(0==Coor[0]%2)&&(Coor[2]==(Coor[0]+1))){
        delete [] Coor;
        return J;
    }
    else{
        delete [] Coor;
        return g;
    }
}

bool check_boundary(int v, int vv){
    if((v%num_leg)<2){
        if((vv/num_leg)>(v/num_leg)){
            return true;
        }
        else{
            return false;
        }
    }
    else{
        if((vv/num_leg)<(v/num_leg)){
            return true;
        }
        else{
            return false;
        }
    }
}

int propagate(string s, int l, vector<int> vStr, int* Spin){
    if(1==(vStr[l]%num_ope)){
        int b=vStr[l]/num_ope;
        int* Coor=new int[num_leg];
        locate_bond(s, b, Coor);
        int s1=search(s, Coor[0], Coor[1]);
        int s2=search(s, Coor[2], Coor[3]);
        Spin[s1]=-Spin[s1];
        Spin[s2]=-Spin[s2];
        delete [] Coor;
    }
    return 1;
}

double mean(double* Dat, int n){
    double s=0.0;
    for(int i=0; i<n; i++){
        s += Dat[i];
    }
    return (s/n);
}

double std_err(double* Dat, int n){
    double s=0.0;
    double a=mean(Dat, n);
    for(int i=0; i<n; i++){
        s += (Dat[i]-a)*(Dat[i]-a);
    }
    return sqrt(s/((n-1)*n));
}

int info(time_t start, time_t end){
    ofstream file_info("task_info");
    file_info<<"lattice: "<<num_dx<<"*"<<num_dy<<endl;
    file_info<<"g="<<g<<endl;
    file_info<<"beta="<<beta<<endl;
    file_info<<"sample: "<<num_sample<<endl;
    file_info<<"time: "<<(end-start)/60.0<<" min"<<endl;
    file_info.close();
    return 1;
}
