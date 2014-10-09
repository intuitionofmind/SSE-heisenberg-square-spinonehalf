/*
 * =====================================================================================
 *
 *       Filename: main.cpp
 *    Description: 
 *        Version: 1.0
 *        Created: 08/22/2014 01:34:53 PM
 *  Last Modified: 10/10/2014 06:33:03 AM
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

int counter = 0;
ofstream file_log("log", ios_base::app | ios_base::binary);

int main(){
    string con = "PBC";
    int ns = num_site(con);
    int n = 0;
    int m = 10;
    time_t start, end;
    int* Spin = new int[ns];
    vector<int> vStr(m, -1);
    start = time(NULL);
    initialize_spin(ns, Spin);
    thermal(con, n, vStr, Spin);
    measure(con, n, vStr, Spin);
    end = time(NULL);
    info(start, end);
    delete [] Spin;
    return 1;
}

bool thermal(string s, int & n, vector<int> & vStr, int* Spin){
    int m = vStr.size();
    for(int i = 0; i<num_thermal; i++){
        diagonal(s, n, vStr, Spin);
        loop(s, n, vStr, Spin);
        int mm = adjust(n, m);
        if(mm>m){
            for(int j = m; j<mm; j++){
                vStr.push_back(-1);
            }
            m = mm;
        }
        file_log<<counter<<endl;
        cout<<counter<<" "<<n<<endl;
        counter++;
    }
    return true;
}

bool measure(string s, int & n, vector<int> & vStr, int* Spin){
    int m = vStr.size();
    int ns = num_site(s);
    int nsc = (num_dy+1)*num_dy;
    ofstream file_n("n.dat", ios_base::app | ios_base::binary);
    double* Ave = new double[ns];
    double* Dat = new double[num_sample*ns];
    double* Tem = new double[num_sample*ns];
    double* pTem = new double[ns*ns];
    double* cTem = new double[nsc*nsc];
    double* Tem_sample = new double[num_sample];
    for(int i = 0; i<num_sample; i++){
        int suc = 0;
        while(suc<num_interval){
            diagonal(s, n, vStr, Spin);
            loop(s, n, vStr, Spin);
            int mm = adjust(n, m);
            if(mm>m){
                for(int j = 0; j<mm; j++){
                    vStr.push_back(-1);
                }
                m = mm;
            }
            suc++;
        }
        for(int j = 0; j<ns; j++){
            Dat[i*ns+j] = 0.5*Spin[j];  //Pauli matrices multiply 1/2 give the Spin operator 
        }
        file_n.write((char*)(&n), sizeof(n));
        file_log<<counter<<endl;
        cout<<counter<<" "<<n<<endl;
        counter++;
    }
    if("PBC" == s){
        ofstream file_c("cor.dat", ios_base::app | ios_base::binary);
        ofstream file_cc("cor_err.dat", ios_base::app | ios_base::binary);
        ofstream file_ne("neel.dat", ios_base::app | ios_base::binary);
        for(int k = 0; k<ns; k++){
            for(int j = 0; j<num_sample; j++){
                for(int i = 0; i<ns; i++){
                    Tem[j*ns+i] = Dat[j*ns+k]*Dat[j*ns+i];
                }
            }
            for(int j = 0; j<ns; j++){
                for(int i = 0; i<num_sample; i++){
                    Tem_sample[i] = Tem[i*ns+j];
                }
                double c = mean(Tem_sample, num_sample);
                double cc = std_err(Tem_sample, num_sample);
                file_c.write((char*)(&c), sizeof(c));
                file_cc.write((char*)(&cc), sizeof(cc));
            }
        }
        for(int k = 0; k<num_sample; k++){
            for(int j = 0; j<ns; j++){
                for(int i = 0; i<ns; i++){
                    pTem[j*ns+i] = Dat[k*ns+j]*Dat[k*ns+i];
                }
            }
            double neel_squa = 0;
            for(int y = 0; y<num_dy; y++){
                for(int x = 0; x<num_dx; x++){
                    int j = search(s, x, y);
                    for(int yy = 0; yy<num_dy; yy++){
                        for(int xx = 0; xx<num_dx; xx++){
                            int i = search(s, xx, yy);
                            neel_squa += cos(PI*(x-xx+y-yy))*pTem[j*ns+i];
                        }
                    }
                }
            }
            neel_squa = 3*neel_squa/double(ns*ns);
            file_ne.write((char*)(&neel_squa), sizeof(neel_squa));
        }
        file_c.close();
        file_cc.close();
        file_ne.close();
    }
    if("CBC" == s){
        ofstream file_ne("neel.dat", ios_base::app | ios_base::binary);
        for(int k = 0; k<num_sample; k++){
            for(int y = 0; y<num_dy; y++){
                for(int x = (num_dx-num_dy)/2; x<(num_dx/2-num_dy/2+num_dy+1); x++){
                    int j = search(s, x, y);
                    int jj = y*(num_dy+1)+(x-(num_dx/2-num_dy/2));
                    for(int yy = 0; yy<num_dy; yy++){
                        for(int xx = (num_dx-num_dy)/2; xx<(num_dx/2-num_dy/2+num_dy+1); xx++){
                            int i = search(s, xx, yy);
                            int ii = yy*(num_dy+1)+(xx-(num_dx/2-num_dy/2));
                            cTem[jj*nsc+ii] = Dat[k*ns+j]*Dat[k*ns+i];
                        }
                    }
                }
            }
            double neel_squa = 0.0;
            for(int y = 0; y<num_dy; y++){
                for(int x = 0; x<(num_dy+1); x++){
                    int j = y*(num_dy+1)+x;
                    for(int yy = 0; yy<num_dy; yy++){
                        for(int xx = 0; xx<(num_dy+1); xx++){
                            int i = yy*(num_dy+1)+xx;
                            neel_squa += cos(PI*(x-xx+y-yy))*cTem[j*nsc+i];
                        }
                    }
                }
            }
            neel_squa = 3*neel_squa/double(nsc*nsc);
            file_ne.write((char*)(&neel_squa), sizeof(neel_squa));
        }
        file_ne.close();
    }
    file_n.close();
    delete [] Ave;
    delete [] Dat;
    delete [] pTem;
    delete [] cTem;
    delete [] Tem_sample;
    return true;
}
