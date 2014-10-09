/*
 * =====================================================================================
 *
 *       Filename: update.cpp
 *    Description: 
 *        Version: 1.0
 *        Created: 08/22/2014 10:48:37 PM
 *  Last Modified: 09/22/2014 11:05:10 PM
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

bool diagonal(string s, int & n, vector<int> & vStr, int* Spin){
    int m = vStr.size();
    int nb;
    if("PBC" == s){
        nb = num_dx*num_dy+num_dy*num_dx;
    }
    else{
        if("CBC" == s){
            nb = num_dx*num_dy+(num_dx+1)*num_dy;
        }
        else{
            return false;
        }
    }
    int* Coor = new int[num_leg];
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis_prob(0, 1);
    uniform_int_distribution<int> dis_bond(0, nb-1);
    for(int l = 0; l<m; l++){
        if(-1 == vStr[l]){
            double r = dis_prob(gen);
            int b = dis_bond(gen);
            locate_bond(s, b, Coor);
            int s1 = search(s, Coor[0], Coor[1]);
            int s2 = search(s, Coor[2], Coor[3]);
            if(Spin[s1] == Spin[s2]){
                continue;
            }
            else{
                if(r<(beta*nb*check_bond(s, b)/(2.0*(m-n)))){
                    vStr[l] = b*num_ope;
                    n++;
                }
            }
        }
        else{
            double r = dis_prob(gen);
            int b = vStr[l]/num_ope;
            int b_type = vStr[l]%num_ope;
            locate_bond(s, b, Coor);
            int s1 = search(s, Coor[0], Coor[1]);
            int s2 = search(s, Coor[2], Coor[3]);
            switch(b_type){
                case 0:
                    if(r<(2.0*(m-n+1)/(beta*nb*check_bond(s, b)))){
                        vStr[l] = -1;
                        n--;
                    }
                    break;
                case 1:
                    Spin[s1] = -Spin[s1];
                    Spin[s2] = -Spin[s2];
                    break;
                default:
                    break;
            }
        }
    }
    delete [] Coor;
    return true;
}


bool off_diagonal(string s, int n, vector<int> & vStr, int* Spin){
    int m = vStr.size();
    int ns = num_site(s);
    vector<int> vStrr(m, 1);
    int* Coor = new int[num_leg];
    int* X = new int[m*num_leg];
    int* First = new int[ns];
    int* Spin_c = new int[ns];
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis_prob(0, 1);
    vertex_list(s, vStr, First, X);
    for(int l = 0; l<m*num_leg; l++){
        if(X[l]<0){
            continue;
        }
        else{
            if(((X[l]^1) == X[l^1])&&(vStr[l/num_leg]%num_ope == vStr[X[l]/num_leg]%num_ope)){
                double r = dis_prob(gen);
                if(r<0.5){
                    vStr[l/num_leg] = (vStr[l/num_leg])^1;
                    vStr[X[l]/num_leg] = (vStr[X[l]/num_leg])^1;
                    X[X[l]] = -2;
                    X[l] = -2;
                    X[X[l^1]] = -2;
                    X[l^1] = -2;
                }
                else{
                    X[X[l]] = -1;
                    X[l] = -1;
                    X[X[l^1]] = -1;
                    X[l^1] = -1;
                }
            }
        }
    }
    /*for(int l = 0; l<m; l++){
        if(-1 == vStr[l]){
            int b = vStr[l]/num_ope;
            int b_type = vStr[l]%num_ope;
            int ll = VerList[l*num_leg+2]/num_leg;
            int lll = VerList[l*num_leg+3]/num_leg;
            if(lll == ll){
                int bb = vStr[ll]/num_ope;
                int bb_type = vStr[ll]%num_ope;
                if((bb_type == b_type)&&(vStrr[l])&&(vStrr[ll])){
                    if(ll>l){
                        vStr[l] = (vStr[l])^1;
                        vStr[ll] = (vStrr[ll])^1;
                        vStrr[l] = 0;
                        vStrr[ll] = 0;
                    }
                    if(ll<l){
                        vStr[l] = (vStr[l])^1;
                        vStr[ll] = (vStr[l])^1;
                        locate_bond(s, b, Coor);
                        int s1 = search(s, Coor[0], Coor[1]);
                        int s2 = search(s, Coor[2], Coor[3]);
                        Spin[s1] = -Spin[s1];
                        Spin[s2] = -Spin[s2];
                        vStrr[l] = 0;
                        vStrr[ll] = 0;
                    }
                }
            }
        }
        else{
            continue;
        }
    }
    for(int i = 0; i<ns; i++){
        Spin_c[i] = Spin[i];
    }
    for(int l = 0; l<m; l++){
        propagate(s, l, vStr, Spin);
    }
    for(int i = 0; i<ns; i++){
        //cout<<(Spin[i]-Spin_c[i])<<" ";
        cout<<(Spin[i])<<" ";
    }
    cout<<endl;*/
    for(int i = 0; i<ns; i++){
        if(-1 == First[i]){
            double r = dis_prob(gen);
            if(r<0.5){
                Spin[i] = -Spin[i];
            }
            continue;
        }
        else{
            if(-2 == X[First[i]]){
                Spin[i] = -Spin[i];
            }
        }
    }
    delete [] Spin_c;
    delete [] Coor;
    delete [] X;
    delete [] First;
    return true;
}

bool loop(string s, int n, vector<int> & vStr, int* Spin){
    int m = vStr.size();
    int ns = num_site(s);
    int* X = new int[m*num_leg];
    int* First = new int[ns];
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis_prob(0, 1);
    vertex_list(s, vStr, First, X);
    for(int l = 0; l<m*num_leg; l++){
        if(X[l]<0){
            continue;
        }
    }
    for(int l = 0; l<m*num_leg; l++){
        if(X[l]<0){
            continue;
        }
        else{
            int ll = l;
            int lll;
            double r = dis_prob(gen);
            if(r<0.5){
                vStr[ll/num_leg] = (vStr[ll/num_leg])^1;
                X[ll] = -2;
                lll = ll^1;
                ll = X[lll];
                X[lll] = -2;
                while(ll != l){
                    vStr[ll/num_leg] = (vStr[ll/num_leg])^1;
                    X[ll] = -2;
                    lll = ll^1;
                    ll = X[lll];
                    X[lll] = -2;
                }
            }
            else{
                X[ll] = -1;
                lll = ll^1;
                ll = X[lll];
                X[lll] = -1;
                while(ll != l){
                    X[ll] = -1;
                    lll = ll^1;
                    ll = X[lll];
                    X[lll] = -1;
                }
            }
        }
    }
    for(int i = 0; i<ns; i++){
        if(-1 == First[i]){
            double r = dis_prob(gen);
            if(r<0.5){
                Spin[i] = -Spin[i];
            }
            continue;
        }
        else{
            if(-2 == X[First[i]]){
                Spin[i] = -Spin[i];
            }
        }
    }
    delete [] X;
    return true;
}

int vertex_list(string s, vector<int> vStr, int* First, int* X){
    int m = vStr.size();
    int ns = num_site(s);
    int* Last = new int[ns];
    int* Coor = new int[num_leg];
    for(int i = 0; i<m*num_leg; i++){
        X[i] = -1;
    }
    for(int i = 0; i<ns; i++){
        First[i] = -1;
        Last[i] = -1;
    }
    for(int l = 0; l<m; l++){
        if(-1 == vStr[l]){
            continue;
        }
        int b = vStr[l]/num_ope;
        locate_bond(s, b, Coor);
        int s1 = search(s, Coor[0], Coor[1]);
        int s2 = search(s, Coor[2], Coor[3]);
        if(s1>s2){
          int s0 = s1;
          s1 = s2;
          s2 = s0;
          }
        int v0 = l*num_leg;
        int v1 = Last[s1];
        int v2 = Last[s2];
        if(-1 != v1){
            X[v1] = v0;
            X[v0] = v1;
        }
        else{
            First[s1] = v0;
        }
        if(-1 != v2){
            X[v2] = v0+1;
            X[v0+1] = v2;
        }
        else{
            First[s2] = v0+1;
        }
        Last[s1] = v0+2;
        Last[s2] = v0+3;
    }
    int f, l;
    for(int i = 0; i<ns; i++){
        f = First[i];
        if(-1 != f){
            l = Last[i];
            X[f] = l;
            X[l] = f;
        }
    }
    delete [] Last;
    delete [] Coor;
    return 1;
}

int adjust(int n, int m){
    int mm = n+n/par_cutoff;
    if(mm<m){
        return m;
    }
    else{
        return mm;
    }
}
