/*
 * =====================================================================================
 *
 *       Filename: head.h
 *    Description: 
 *        Version: 1.0
 *        Created: 08/22/2014 01:32:22 PM
 *  Last Modified: 09/21/2014 08:39:37 PM
 *       Compiler: gcc
 *
 *         Author: Wei Zheng
 *          Email: intuitionofmind@gmail.com
 *   Organization: Institute for Advanced Study, Tsinghua University
 *
 * =====================================================================================
 */


#include <iostream>
#include <fstream>
#include <random>
#include<iomanip>
#include <memory.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <malloc.h>
#include <ctime>
#include <algorithm>
#include <initializer_list>
#include <vector>
#include <map>
//#include <mpi.h>

#define PI 3.14159265358979

using namespace std;

extern int counter;
int main();
bool thermal(string s, int & n, vector<int> & vStr, int* Spin);
bool measure(string s, int & n, vector<int> & vStr, int* Spin);

int num_site(string s);
double unit_prob();
int search(string c, int x, int y);
int initialize_spin(int len, int* Spin);
bool locate_bond(string s, int b, int* Coor);
double check_bond(string s, int b);
bool check_boundary(int v, int vv);
int propagate(string s, int l, vector<int> vStr, int* Spin);
double mean(double* Dat, int n);
double std_err(double* Dat, int n);
int info(time_t start, time_t end);

bool diagonal(string s, int & n, vector<int> & vStr, int* Spin);
bool off_diagonal(string s, int n, vector<int> & vStr, int* Spin);
bool off_diagonal2(string s, int n, vector<int> & vStr, int* Spin);
bool loop(string s, int n, vector<int> & vStr, int* Spin);
int vertex_list(string s, vector<int> vStr, int* First, int* X);
int adjust(int n, int m);
