/*
 * =====================================================================================
 *
 *       Filename: main.cpp
 *    Description: 
 *        Version: 1.0
 *        Created: 08/22/2014 01:34:53 PM
 *  Last Modified: 02/25/2015 10:16:45 AM
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

int counter=0;
ofstream file_log("log", ios_base::app | ios_base::binary);

int main(){
    string con="CBC";
    int ns=num_site(con);
    int n=0;
    int m=10;
    time_t start, end;
    int* Spin=new int[ns];
    vector<int> vStr(m, -1);
    start=time(NULL);
    initialize_spin(ns, Spin);
    thermal(con, n, vStr, Spin);
    measure(con, n, vStr, Spin);
    end=time(NULL);
    info(start, end);
    delete [] Spin;
    return 1;
}

bool thermal(string s, int & n, vector<int> & vStr, int* Spin){
    int m=vStr.size();
    for(int i=0; i<num_thermal; i++){
        diagonal(s, n, vStr, Spin);
        loop(s, n, vStr, Spin);
        int mm=adjust(n, m);
        if(mm>m){
            for(int j=m; j<mm; j++){
                vStr.push_back(-1);
            }
            m=mm;
        }
        file_log<<counter<<endl;
        counter++;
    }
    return true;
}

bool measure(string s, int & n, vector<int> & vStr, int* Spin){
    int m=vStr.size();
    int ns=num_site(s);
    double* sCor=new double[ns*ns];
    int* Coor=new int[num_leg];
    ofstream file_ne("neel.dat", ios_base::app | ios_base::binary);
    for(int i=0; i<num_sample; i++){
        int suc=0;
        while(suc<num_interval){
            diagonal(s, n, vStr, Spin);
            loop(s, n, vStr, Spin);
            int mm=adjust(n, m);
            if(mm>m){
                for(int j=0; j<mm; j++){
                    vStr.push_back(-1);
                }
                m=mm;
            }
            suc++;
        }
        for(int j=0; j<ns; j++){
            for(int k=0; k<ns; k++){
                sCor[j*ns+k]=(0.5*Spin[j])*(0.5*Spin[k]);
            }
        }
        double neel_squa=0.0;
        if("PBC"==s){
            for(int y=0; y<num_dy; y++){
                for(int x=0; x<num_dx; x++){
                    int j=search(s, x, y);
                    for(int yy=0; yy<num_dy; yy++){
                        for(int xx=0; xx<num_dx; xx++){
                            int i=search(s, xx, yy);
                            neel_squa+=cos(PI*(x-xx+y-yy))*sCor[j*ns+i];
                        }
                    }
                }
            }
            neel_squa=3*neel_squa/double(ns*ns);
        }
        if("CBC"==s){
            if(0==((num_dx-num_dy)/2)%2){  //note that unit cell along x direction has been doubled and may cause something wrong
                for(int y=0; y<num_dy; y++){
                    for(int x=(num_dx-num_dy)/2; x<(num_dx/2-num_dy/2+num_dy); x++){
                        int j=search(s, x, y);
                        for(int yy=0; yy<num_dy; yy++){
                            for(int xx=(num_dx-num_dy)/2; xx<(num_dx/2-num_dy/2+num_dy); xx++){
                                int i=search(s, xx, yy);
                                neel_squa+=cos(PI*(x-xx+y-yy))*sCor[j*ns+i];
                            }
                        }
                    }
                }
                int nsc=num_dy*num_dy;
                neel_squa=3*neel_squa/double(nsc*nsc);
            }
            else{
                for(int y=0; y<num_dy; y++){
                    for(int x=(num_dx/2-num_dy/2+1); x<(num_dx/2-num_dy/2+num_dy+1); x++){
                        int j=search(s, x, y);
                        for(int yy=0; yy<num_dy; yy++){
                            for(int xx=(num_dx-num_dy)/2; xx<(num_dx/2-num_dy/2+num_dy); xx++){
                                int i=search(s, xx, yy);
                                neel_squa+=cos(PI*(x-xx+y-yy))*sCor[j*ns+i];
                            }
                        }
                    }
                }
                int nsc=num_dy*num_dy;
                neel_squa=3*neel_squa/double(nsc*nsc);
            }
        }
        file_ne.write((char*)(&neel_squa), sizeof(neel_squa));
        file_log<<counter<<endl;
        counter++;
    }
    delete [] sCor;
    delete [] Coor;
    file_ne.close();
    return true;
}
