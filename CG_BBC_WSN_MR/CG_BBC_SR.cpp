//g++ -m64 -g -o CG_BBC_SR.out CG_BBC_SR.cpp -I/opt/gurobi562/linux64/include -L/opt/gurobi562/linux64/lib -lgurobi_c++ -L/opt/gurobi562/linux64/lib -lgurobi56 -lpthread -lm

#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <vector>
#include <set>
#include "gurobi_c++.h"
#include <fstream>
#include <math.h>
#include <limits.h>
using namespace std;
 
 
//Data structure definitions
typedef vector<int> Vec;
typedef vector<Vec> Mat;
typedef vector<Mat> Array3D;
 
 
//Declarations
double metodos(string dir,   int method,double Rs,double Rc);
double Column_Generation_ILP(string dir,  int method,double Rs,double Rc);
double Column_Generation_CP(string dir,  int method,double Rs,double Rc);
 
int ConnBin(double x1, double x2, double y1, double y2,double Rs);
int Check_cov(Mat &Adj, Vec const &Solution, int T, int S); 
 
 
//Benders declarations
vector<double> Benders_Flow(Mat &A1, Mat &y, int Nodes, int Targets, double *Bound,int *status);
vector<double> Benders_Algorithm(int Nodes, int Targets,int Nroles,  double dual[], Mat &A1,vector<double> &nrj,double *ObjAux);
vector<double> Benders_Algorithm_Complete(int Nodes, int Targets,int Nroles,  double dual[], Mat &A1,vector<double> &nrj,double *ObjAux);
 
//Functions used to create the models
GRBModel Gen_Master(int S);
GRBModel Gen_Aux_MultiRole(int Nodes, int Targets,int Nroles, Mat &A1,vector<double> &nrj);
GRBModel Gen_Cutting_planes(int Nodes, int Targets,int Nroles,   Mat const &A1,vector<double> &nrj);
 
//Functions used to identify connected components
Mat Con_Comp(Mat &Adj, Vec &Solution);
void DFS(Mat &Adj, Vec const &Solution, Vec &labels, int cur_node, Vec &Comp);
 
//Callback class definition
vector<double> Cutting_planes(int Nodes, int Targets,int Nroles,  double dual[], Mat &A1,vector<double> &nrj,double *ObjAux);
 
 
//Callback related stuff
vector< vector<int> > NewCols;
//global model used to manage BSP
GRBEnv* envBenders = new GRBEnv(); //creation of model environment
GRBModel BSP = GRBModel(*envBenders);
 
 
class connectivity: public GRBCallback
{
public:
    GRBVar* vars;
    int n_nodes;
    int nrj_lev;
    int targets;
    vector< vector<int> > Adj;
    vector< vector<int> > Cols;
    connectivity(GRBVar* Y_vars, int xn, int T,vector< vector<int> > const &A1) {
        vars = Y_vars;
        n_nodes = xn;
        Adj=A1;
        targets=T;
    }
protected:
    void callback() {
        try {
            if (where == GRB_CB_MIPSOL) {
                 
                // Found an integer feasible solution - does it visit every node?
                vector<int> x;
                vector<int> members;
                int coverage=0;
                Mat CC;
                for (int i = 0; i < n_nodes; i++){
                    //The numbers one and two are used to denote the index of the Y solution according to the representation in the solutions
                    x.push_back(round(getSolution(vars[i])));
                }
                //x[0]=1;
                
                CC=Con_Comp(Adj, x);
                for (int i=0; i<CC[0].size(); i++) {
					if (x[CC[0][i]]>0) {
						members.push_back(CC[0][i]);
					}
				}
				coverage=Check_cov(Adj, members,targets, n_nodes-1);
				if (coverage==targets) {
					Vec SaveMe;
					for (int i = 0; i < n_nodes; i++){
						//The numbers one and two are used to denote the index of the Y solution according to the representation in the solutions
						SaveMe.push_back(0);
					}
					for (int i = 0; i < CC[0].size(); i++){
						//The numbers one and two are used to denote the index of the Y solution according to the representation in the solutions
						SaveMe[CC[0][i]]=x[CC[0][i]];
					}
					NewCols.push_back(SaveMe);
				}

                if (CC.size()>1 and coverage<targets) {
                    // Add connectiivity elimination constraint
                    for (int i=1; i<CC.size(); i++) {
                        GRBLinExpr left = 0;
                        GRBLinExpr right =0;
                        for (int j=0; j<CC[i].size(); j++) {
                            left+=vars[CC[i][j]];
                            for (int k=0; k<n_nodes; k++) {
                                if ((Adj[CC[i][j]][k]>0 or Adj[k][CC[i][j]]>0) and (x[k]==0)) {
                                    right+=vars[k];
                                }
                            }
                        }
                        addLazy(left<=CC[i].size()*right);
                        //cout<<left<<"<="<<CC[i].size()<<"*"<<right<<endl;
                         
                    }
                    //Combinatorial cut
                    //*
                    GRBLinExpr left = 0;
                    for (int i=1; i<n_nodes; i++) {
                        if ((x[i]==0)) {
                            left+=vars[i];
                        }
                        else{
                            left+=1-vars[i];
                             
                        }
                    }
                    addLazy(left>=1);//*/
                }
            }
        } catch (GRBException e) {
            cout << "Error number: " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) {
            cout << "Error during callback" << endl;
        }
    }
};
 
 
//Main program
int main(){
    srand (time(0));
    string sensores="arbol", targets="_verde";
    int nsens[]={100,200,300,400,500};
    int ntar[]={15,30};
    int nBS[]={1};
    double RadS[]={100,125};
    double RadC[]={125,125};
     
    string dir;
    double inicio,final;
    double MLP;
    /*
     for (int m=1; m<=1; m++) {
     for (int size=0; size<=4; size++) {
     
     for (int Stations=0; Stations<=0; Stations++) {
     
	 for (int targets=0; targets<=1; targets++) {
     for (int Instance=0; Instance<=3; Instance++) {
     stringstream S,T,BS,ins;
     dir="/home/fabian/Documents/MLP/Instance_";
     S<<nsens[size];
     T<<ntar[targets];
     BS<<nBS[Stations];
     ins<<Instance;
     dir+=BS.str();
     dir+="_";
     dir+=S.str();
     dir+="_";
     dir+=T.str();
     dir+="_";
     dir+=ins.str();
     dir+=".txt";
     cout<<dir<<endl;
     //inicio=clock();
     MLP=Column_Generation_CP(dir,3,RadS[m],RadC[m]);
     //final=clock();
      
      
     }
     }
     }
     }
     }//*/
     
    time_t start,end;
    time (&start);
    dir="/home/fabian/Documents/MLP/Instance_1_200_15_1.txt";
    MLP=Column_Generation_CP(dir,2,125,125);
    //MLP=Stab_CG_CP(dir,2,100,125);
    //MLP=Neame_Generation_CP(dir,2,100,125);
    //MLP=Column_Generation_ILP(dir,2,100,125);
     
    time (&end);
    cout<<" time: "<< difftime (end,start) <<endl;//*/
     
}
 
 
double Column_Generation_CP(string dir,  int method,double Rs,double Rc){
    ofstream FILE_Exact_Results;
      int cont_sense=0;
    //===========================Starting Data reading==================================================
    ifstream indata; // indata is like cin
    vector<double> coorx;
    vector<double> coory;
    double num; // variable for input value
    indata.open(dir.c_str()); // opens the file
    if(!indata) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    int   S, T, BS;
    indata>> BS;
    indata >> S; // Number of sensors
    indata >> T; // Number of targets
     
    int count=0;
    while ( count<1+S+T ) { // keep reading until end-of-file
        indata >> num; // sets EOF flag if no value found
        coorx.push_back(num);
        indata >> num; // sets EOF flag if no value found
        coory.push_back(num);
        count++;
    }
    indata.close();
     
    //===========================End Data reading========================================================
     
    //Declarations of data input for model and CG approach
    Vec fila;
    Mat A1;
     
     
    //Preparing the matrix related with the edges
     
    for(int i=0; i<1+S+T;i++){
        for(int j=0; j<1+S+T;j++){
            fila.push_back(0);
        }
        A1.push_back(fila);
        fila.clear();
    }
     
     
    //Preparing the matrix related with the connectivity placing 1's where the connectivity requirement exists
    for(int i=0; i<1+S+T;i++){
        for(int j=0; j<1+S+T;j++){
            if(i==0 or i==j){
                A1[i][j]=0;
            }
            if(i>0 and i<1+S){
                if (j<1+S) {
                    A1[i][j]=ConnBin(coorx[i], coorx[j], coory[i], coory[j],Rc);
                }
            }
            if(i>=1+S){
                if(j>0 and j<1+S){
                    A1[i][j]=ConnBin(coorx[i], coorx[j], coory[i], coory[j],Rs);
                }
            }
             
        }
    }
     
 
     
     
    //Preparing the creation of the Model
    int ncols=0; //Variable to store the number of generated columns
    GRBModel MP=Gen_Master(S);
     
    //=================================================================
    GRBEnv* envPS = 0;
    envPS = new GRBEnv(); //creation of model environment
    GRBVar* Y = 0; //variables concerning to the use of the sensors and the related state
     
    GRBModel PS = GRBModel(*envPS); //Declaration ofmodel
    PS.set(GRB_StringAttr_ModelName, "Benders_Master"); //Model Name Declaration
    PS.getEnv().set(GRB_IntParam_OutputFlag, 0);
    PS.getEnv().set(GRB_IntParam_LazyConstraints, 1);
    //PS.getEnv().set(GRB_DoubleParam_Heuristics,0.0);
    PS.getEnv().set(GRB_IntParam_MIPFocus,1);
     
    // The objective is to maximize the total fixed and variable costs
    PS.set(GRB_IntAttr_ModelSense, 1);
    PS.getEnv().set(GRB_DoubleParam_TimeLimit, 10);
    //model.getEnv().set(GRB_IntParam_Presolve,0);
     
     
 
    Y = PS.addVars(S+1,GRB_BINARY);
    PS.update();
    for (int p = 0; p < S+1; ++p)
    {
        ostringstream vname;
        vname << "Y_" << p;
        Y[p].set(GRB_DoubleAttr_Obj, 0);
        Y[p].set(GRB_StringAttr_VarName, vname.str());
    }
     
     
    PS.update();
    //Adding initial cuts
    // Role Allocation to Sensors
 
    for (int v = S+1;v < S+1+T; v++)
    {
        GRBLinExpr User_Cuts = 0;
        for (int u = 1; u < S+1; u++){
            User_Cuts += Y[u]*A1[v][u];
        }
        ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
        cname << "User_Cuts_" << v;
        PS.addConstr(User_Cuts >= 1, cname.str());
    }
     
    for (int v = 0;v < 1; v++)
    {
        GRBLinExpr User_Cuts = 0;
        for (int u = 1; u < S+1; u++){
            User_Cuts += (Y[u])*A1[u][v];
        }
        ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
        cname << "User_Cuts_" << v;
        PS.addConstr(User_Cuts >= 1, cname.str());
    }
     
    for (int v = 1;v < S+1; v++)
    {
        if(A1[v][0]<=0){
            GRBLinExpr User_Cuts = 0;
            for (int u = 1; u < S+1; u++){
                User_Cuts += (Y[u])*A1[u][v];
            }
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "User_Cuts_" << v;
            PS.addConstr(User_Cuts >= Y[v], cname.str());
        }
         
    }
     
    GRBLinExpr User_Cuts = 0;
    for (int v = 1;v < S+1; v++)
    {
        User_Cuts += (Y[v]);
    }
    ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
    cname << "Max_relays";
    int RHS=25,Remo_ctr=0;
    PS.addConstr(User_Cuts <= RHS, cname.str());
     
    PS.update();
    //The model is solved calling the "callback" used to add the connectivity constraints
    //model.write("/Users/Fabian/Dropbox/Imprimid.lp");
     
    //solution pool: is used to sttore the solutions found along the CG process
    connectivity cb=connectivity(Y, S+1,T, A1);
    PS.setCallback(&cb);
     
    //Decaring constraints
    GRBConstr* c = 0;
    GRBConstr* Ps_ctr = 0;
     
    c=MP.getConstrs();
     
    Ps_ctr=PS.getConstrs();
    //Initializing the number of columns
    ncols=1;
     
    time_t start,end;
    time (&start);
    int iters=0;
    
    //initialization
  
	//PS.update();
	//PS.write("/Users/Fabian/Desktop/Salida.lp");
	NewCols.clear();
	PS.getEnv().set(GRB_DoubleParam_Cutoff,1);

    while (true) {
        iters+=1;
        MP.optimize();
        time(&end);
        cout<<" Lifetime:"<<MP.get(GRB_DoubleAttr_ObjVal)<<" Time: "<<difftime(end, start)<<" n cols: " << ncols<<" n iters: " << iters<<endl;
         
        //This is a subroutine used to modify the objective function of the PS according to the duals found throught the MP!
        for(int i=0;i<S;i++){
                ostringstream sname;
                //The constraint i is related to the variables of the Nodes i+1
                sname << "Y_" << i+1;
                Y[i+1].set(GRB_DoubleAttr_Obj,c[i].get(GRB_DoubleAttr_Pi));
        }
        //PS.update();
        //PS.write("/Users/Fabian/Desktop/Salida.lp");
        NewCols.clear();
        PS.optimize();
        time(&end);
        if(PS.get(GRB_IntAttr_Status)!=3 and 1-PS.get(GRB_DoubleAttr_ObjVal)>0.0000001 and  (PS.get(GRB_IntAttr_Status)==2 or PS.get(GRB_IntAttr_Status)==13) and difftime(end, start)<3600){
            for (int p=NewCols.size()-1; p<NewCols.size(); p++) {
                ncols++;
                //cout <<"Adding Col "<<ncols<<" Lifetime:"<<MP.get(GRB_DoubleAttr_ObjVal)<<endl;
                GRBColumn Col;
                //The columns, in the RMP, are created with the energy of the sensors exactly as they are stacked in the E_x vector, whic starts with the first sensor
                double E_tot=0;
                for (int i=0; i<S; i++) {
                    //The constraint i is related to the variables of the Nodes i+1
                    E_tot=NewCols[p][i+1];
                    //cout<<E_tot<<" ";
                    if (E_tot>0)
                        Col.addTerm(E_tot, c[i]);
                }
                 
                //cout<<endl;
                MP.addVar(0.0, GRB_INFINITY, 1, GRB_CONTINUOUS,Col,"t"+ncols);
                 
            }
            PS.getEnv().set(GRB_DoubleParam_TimeLimit, 10);
            cont_sense=0;
        }
        else{
            if ((PS.get(GRB_IntAttr_Status)!=GRB_OPTIMAL or PS.get(GRB_IntAttr_Status)==3) and Remo_ctr==0 and difftime(end, start)<3600) {
                if (Remo_ctr==0) {
                    cout<<"Modifying time limit "<<endl;
                    PS.getEnv().set(GRB_DoubleParam_TimeLimit, 10*(cont_sense+1));

                    cont_sense+=1;
                }
                if(NewCols.size()>0){
                for (int p=NewCols.size()-1; p<NewCols.size(); p++) {
                    ncols++;
                    //cout <<"Adding Col "<<ncols<<" Lifetime:"<<MP.get(GRB_DoubleAttr_ObjVal)<<endl;
                    GRBColumn Col;
                    //The columns, in the RMP, are created with the energy of the sensors exactly as they are stacked in the E_x vector, whic starts with the first sensor
                    double E_tot=0;
                    for (int i=0; i<S; i++) {
                        //The constraint i is related to the variables of the Nodes i+1
                        E_tot=NewCols[p][i+1];
                        //cout<<E_tot<<" ";
                        if (E_tot>0)
                            Col.addTerm(E_tot, c[i]);
                    }
                     
                    //cout<<endl;
                    MP.addVar(0.0, GRB_INFINITY, 1, GRB_CONTINUOUS,Col,"t"+ncols);
                     
                }}
				if (cont_sense>100) {
                    cout<<"removing extra ctr"<<endl;
                    PS.remove(Ps_ctr[PS.get(GRB_IntAttr_NumConstrs)-1]);
                    Remo_ctr=1;
                }
 
            }
            else{
                    FILE_Exact_Results.open("/home/fabian/Documents/CG_BBC_WSN_MR/Res_CG_BBC_SR_0.txt",fstream::app);
                    FILE_Exact_Results<<dir<<" Lifetime:"<<MP.get(GRB_DoubleAttr_ObjVal)<<" Time: "<<difftime(end, start)<<" n cols: " << ncols<<" n iters: " << iters<<endl;
                    FILE_Exact_Results.close();
                    break;
                
            }
             
        }
         
    }
    return 0;
}
 
 
double Column_Generation_ILP(string dir,  int method,double Rs,double Rc){
    ofstream FILE_Exact_Results;
    //===========================Starting Data reading==================================================
    ifstream indata; // indata is like cin
    vector<double> coorx;
    vector<double> coory;
    double num; // variable for input value
    indata.open(dir.c_str()); // opens the file
    if(!indata) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    int   S, T, BS;
    indata>> BS;
    indata >> S; // Number of sensors
    indata >> T; // Number of targets
     
    int count=0;
    while ( count<1+S+T ) { // keep reading until end-of-file
        indata >> num; // sets EOF flag if no value found
        coorx.push_back(num);
        indata >> num; // sets EOF flag if no value found
        coory.push_back(num);
        count++;
    }
    indata.close();
     
    //===========================End Data reading========================================================
     
    //Declarations of data input for model and CG approach
    Vec fila;
    Mat A1;
    double inicia=clock();
     
    //Preparing the matrix related with the edges
    for(int i=0; i<1+S+T;i++){
        for(int j=0; j<1+S+T;j++){
            fila.push_back(0);
        }
        A1.push_back(fila);
        fila.clear();
    }
     
    //Preparing the matrix related with the connectivity placing 1's where the connectivity requirement exists
    for(int i=0; i<1+S+T;i++){
        for(int j=0; j<1+S+T;j++){
            if(i==0 or i==j){
                A1[i][j]=0;
            }
            if(i>0 and i<1+S){
                if (j<1+S) {
                    A1[i][j]=ConnBin(coorx[i], coorx[j], coory[i], coory[j],Rc);
                }
            }
            if(i>=1+S){
                if(j>0 and j<1+S){
                    A1[i][j]=ConnBin(coorx[i], coorx[j], coory[i], coory[j],Rs);
                }
            }
             
        }
    }
     
     
    vector<double> nrj;
    nrj.push_back(0);
    nrj.push_back(1.0);
    nrj.push_back(1.0);
     
    //Preparing the creation of the Model
    int ncols=0; //Variable to store the number of generated columns
    GRBModel MP=Gen_Master(S);
    GRBModel PS=Gen_Aux_MultiRole(S+1, T, nrj.size(), A1, nrj);
     
    PS.getEnv().set(GRB_IntParam_OutputFlag, 0);
    PS.getEnv().set(GRB_IntParam_LazyConstraints, 1);
    PS.getEnv().set(GRB_DoubleParam_Heuristics,0.0);
    // The objective is to maximize the total fixed and variable costs
    PS.set(GRB_IntAttr_ModelSense, 1);
    PS.getEnv().set(GRB_DoubleParam_TimeLimit, 10);
    //Decaring constraints
    GRBConstr* c = 0;
    c=MP.getConstrs();
     
    GRBConstr* Ps_ctr = 0;
    Ps_ctr=PS.getConstrs();
    //Initializing the number of columns
    ncols=1;
    int Remo_ctr=0,RHS=11;
    time_t start,end;
    time (&start);
    while (true) {
        MP.optimize();
        cout<<" Lifetime:"<<MP.get(GRB_DoubleAttr_ObjVal)<<endl;
         
        //This is a subroutine used to modify the objective function of the PS according to the duals found throught the MP!
        for(int i=0;i<S;i++){
            for (int j=1; j<nrj.size(); j++) {
                ostringstream sname;
                //The constraint i is related to the variables of the Nodes i+1
                sname << "Y_" << i+1 <<"-"<<j;
                PS.getVarByName(sname.str()).set(GRB_DoubleAttr_Obj,c[i].get(GRB_DoubleAttr_Pi));
            }
        }
        //PS.update();
        //PS.write("/Users/Fabian/Desktop/Salida.lp");
        PS.optimize();
        time(&end);
        if(PS.get(GRB_IntAttr_Status)!=3 and 1-PS.get(GRB_DoubleAttr_ObjVal)>0.0000001 and  (PS.get(GRB_IntAttr_Status)==2 or PS.get(GRB_IntAttr_Status)==13)){
            ncols++;
            //cout <<"Adding Col "<<ncols<<" Lifetime:"<<MP.get(GRB_DoubleAttr_ObjVal)<<endl;
            GRBColumn Col;
            //The columns, in the RMP, are created with the energy of the sensors exactly as they are stacked in the E_x vector, whic starts with the first sensor
            for (int i=0; i<S; i++) {
                double E_tot=0;
                for (int j=0; j<3; j++) {
                    ostringstream sname;
                    //The constraint i is related to the variables of the Nodes i+1
                    sname << "Y_" << i+1 <<"-"<<j;
                    E_tot+=PS.getVarByName(sname.str()).get(GRB_DoubleAttr_X)*nrj[j];
                }
                if (E_tot>0) {
                    Col.addTerm(E_tot, c[i]);
                }
            }
            //cout<<"ya"<<endl;
            MP.addVar(0.0, GRB_INFINITY, 1, GRB_CONTINUOUS,Col,"t"+ncols);
            //MP.update();
        }
        else{
            if ((PS.get(GRB_IntAttr_Status)!=GRB_OPTIMAL or PS.get(GRB_IntAttr_Status)==3) and Remo_ctr==0 and difftime(end, start)<3600) {
                if (Remo_ctr==0) {
                    cout<<"Modifying extra ctr"<<endl;
                    PS.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);
                    RHS=RHS*2;
                    Ps_ctr[PS.get(GRB_IntAttr_NumConstrs)-1].set(GRB_DoubleAttr_RHS, RHS);
                }
            }
            else{
                if (PS.get(GRB_IntAttr_Status)==GRB_OPTIMAL and Remo_ctr==0 and 1-PS.get(GRB_DoubleAttr_ObjVal)<=0.0000001 and difftime(end, start)<3600) {
                    cout<<"removing extra ctr"<<endl;
                    PS.remove(Ps_ctr[PS.get(GRB_IntAttr_NumConstrs)-1]);
                    Remo_ctr=1;
                }
                else{
                    FILE_Exact_Results.open("/home/fabian/Documents/CG_BBC_WSN_MR/Res_CG_BBC_SR.txt",fstream::app);
                    FILE_Exact_Results<<dir<<" Lifetime:"<<MP.get(GRB_DoubleAttr_ObjVal)<<" Time: "<<difftime(end, start)<<" n cols: " << ncols<<endl;
                    FILE_Exact_Results.close();
                    break;
                }
            }
             
        }
    }
    return 0;
}
 
int ConnBin(double x1, double x2, double y1, double y2,double Rs){
    double dist;
    double a=0;
     
    dist=sqrt(pow(x1-x2,2)+pow(y1-y2,2));
    if(dist<=Rs){
        a=1;
    }
    return a;
}
 
GRBModel Gen_Aux_MultiRole(int Nodes, int Targets,int Nroles, Mat &A1,vector<double> &nrj){
    //Declaration of environments and variables for the model
    GRBEnv* env = 0;
    GRBVar** y = 0;
    GRBVar** x = 0;
    vector<double> Column;
    env = new GRBEnv(); //creation of model environment
    GRBModel model = GRBModel(*env); //Declaration ofmodel
     
    int count=0;
     
    double FObj;
    try{
         
        // Model
         
        //model.getEnv().set(GRB_DoubleParam_IntFeasTol,1e-9);
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        model.set(GRB_StringAttr_ModelName, "Pricing Subproblem"); //Model Name Declaration
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        //model.getEnv().set(GRB_DoubleParam_Heuristics ,0.0);
        //************************** ******* *************************
         
        //***********************************************************
        //Creating Activation mode variables
        y = new GRBVar* [Nodes];
        for (int w = 0; w < Nodes; ++w) //Three is for the number of activation modes
        {
            y[w] = model.addVars(nrj.size(),GRB_BINARY);
            model.update();
            for (int p = 0; p < nrj.size(); ++p)
            {
                ostringstream vname;
                vname << "Y_" << w << "-" << p;
                y[w][p].set(GRB_StringAttr_VarName, vname.str());
            }
        }
         
        //************************** ******* *************************
        //Creating flow variables
         
        x = new GRBVar* [Nodes+Targets];
        for (int w = 0; w < Nodes+Targets; ++w)
        {
             
            x[w] = model.addVars(Nodes+Targets,GRB_INTEGER);
            count++;
            model.update();
            for (int p = 0; p < Nodes+Targets; ++p)
            {
                ostringstream vname;
                vname << "X_" << w << "_" << p;
                //x[w][p].set(GRB_DoubleAttr_Obj, 0);
                x[w][p].set(GRB_StringAttr_VarName, vname.str());
            }
        }
         
        model.update(); //Updating model to include variables
         
        //*************************************************** Including Balance Constraints  *************************************************
         
        // Flow constraints concerning to nodes
         
        for (int v = 0; v < Nodes+Targets; v++)
        {
            GRBLinExpr Act_Mode = 0;
            if (v==0) {
                for (int u = 0; u < Nodes+Targets; u++)
                {
                    Act_Mode += x[v][u]*A1[v][u];
                    Act_Mode += x[u][v]*A1[u][v]*(-1);
                     
                }
                ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
                cname << "Balance_" << v;
                model.addConstr(Act_Mode == -Targets, cname.str());
            }
            else{if (v>0 and v<Nodes) {
                for (int u = 0; u < Nodes+Targets; u++)
                {
                    Act_Mode += x[v][u]*A1[v][u];
                    Act_Mode += x[u][v]*A1[u][v]*(-1);
                     
                }
                ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
                cname << "Balance_" << v;
                model.addConstr(Act_Mode == 0, cname.str());
            }
            else{
                 
                for (int u = 0; u < Nodes+Targets; u++)
                {
                    Act_Mode += x[v][u]*A1[v][u];
                    Act_Mode += x[u][v]*A1[u][v]*(-1);
                     
                }
                ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
                cname << "Balance" << v;
                model.addConstr(Act_Mode == 1, cname.str());
                 
            }
            }
             
        }
         
        //*************************************************** Including Energy Level Constraint  *************************************************
         
        // Role Allocation to Sensors
        for (int v = 1;v < Nodes; v++)
        {
            GRBLinExpr Role_Alloc = 0;
            for (int l = 0; l < 3; l++){
                Role_Alloc += y[v][l];
            }
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "Role_Alloc_" << v;
            model.addConstr(Role_Alloc == 1, cname.str());
        }
         
         
        for (int v = Nodes;v < Nodes+Targets; v++)
        {
            GRBLinExpr User_Cuts = 0;
            for (int u = 1; u < Nodes; u++){
                User_Cuts += y[u][2]*A1[v][u];
            }
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "User_Cuts_" << v;
            model.addConstr(User_Cuts >= 1, cname.str());
        }
         
        for (int v = 0;v < 1; v++)
        {
            GRBLinExpr User_Cuts = 0;
            for (int u = 1; u < Nodes; u++){
                User_Cuts += (y[u][2]+y[u][1])*A1[u][v];
            }
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "User_Cuts_" << v;
            model.addConstr(User_Cuts >= 1, cname.str());
        }
        //*************************************************** Coverage Constraints  *************************************************
         
        //Including coverage Constraints
        for (int u = Nodes;u < Nodes+Targets; u++)
        {
            GRBLinExpr Cov_Constraint = 0;
            for (int v = 1; v < Nodes; v++){
                Cov_Constraint += y[v][2]*A1[u][v];
            }
             
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "Cov_Constraint" << u;
            model.addConstr(Cov_Constraint >= 1, cname.str());
        }
         
         
        //*************************************************** Flow Bounds based on connectivit and coverage  *****************************************
         
        // Role Allocation to Sensors
        for (int v = 0;v < Nodes; v++)
        {
            for (int u = 0; u < Nodes; u++){
                GRBLinExpr FlowBounds = 0;
                for (int l=1; l<3; l++) {
                    FlowBounds += y[v][l];
                }
                ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
                cname << "Bound_Sensor_sensor" << u<<"_"<<v;
                model.addConstr(A1[u][v]*x[u][v] <= FlowBounds*Targets, cname.str());
            }
        }
         
        for (int v = 0;v < Nodes; v++)
        {
            for (int u = 0; u < Nodes; u++){
                GRBLinExpr FlowBounds = 0;
                for (int l=1; l<3; l++) {
                    FlowBounds += y[v][l];
                }
                ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
                cname << "Bound_Sensor_sensor" << u<<"_"<<v;
                model.addConstr(A1[v][u]*x[v][u] <= FlowBounds*Targets, cname.str());
            }
        }
         
        for (int v = 0;v < Nodes; v++)
        {
            for (int u = Nodes; u < Nodes+Targets; u++){
                GRBLinExpr FlowBounds = 0;
                FlowBounds += y[v][2];
                ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
                cname << "Bound_Target_sensor" << u<<"_"<<v;
                model.addConstr(A1[u][v]*x[u][v] <= FlowBounds, cname.str());
            }
        }
         
         
        //Fake constraint to stabilize
        int RHS=10;
        GRBLinExpr User_Cuts = 0;
        for (int v = 1;v < Nodes; v++)
        {
            User_Cuts += (y[v][1]+y[v][2]);
        }
        ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
        cname << "Max_relays";
        model.addConstr(User_Cuts <= RHS, cname.str());
         
        model.update();
    }
    catch(...){
        cout<<"Error in Auxiliary problem creation";
         
    }
    return model;
     
}
 
GRBModel Gen_Master(int S){
    //Environment and model related variables
    GRBEnv* envMP = 0;
    GRBVar* t_var = 0;
     
    //Initialization
    vector< vector<double> > Sets;
    vector<double> ColIni;
    for (int i=0; i<S; i++) {
        ColIni.push_back(1);
    }
     
    Sets.push_back(ColIni);
    ColIni.clear();
     
    // Model
    envMP = new GRBEnv();
    GRBModel master = GRBModel(*envMP);
    master.set(GRB_StringAttr_ModelName, "MR_CMLP");
    // Plant open decision variables: open[p] == 1 if plant p is open.
    t_var = master.addVars(S, GRB_CONTINUOUS);
    master.update();
    master.getEnv().set(GRB_IntParam_OutputFlag, 0);
     
    //Decalaration of model's variables for the master problem
    for (int w = 0; w < 1; w++)
    {
        ostringstream vname;
        vname << "t" << w;
        t_var[w].set(GRB_DoubleAttr_Obj, 1); //el 1 indica que solo cuesta 1 en la fuciÃ³n objetivo usarla
        t_var[w].set(GRB_StringAttr_VarName, vname.str());
         
    }
     
    // The objective is to maximize the total fixed and variable costs
    master.set(GRB_IntAttr_ModelSense, -1);
     
    // Update model to integrate new variables
    master.update();
     
    // Battery constraints
    // A new constraint for each sensor is added
     
    for (int p = 0; p < S; p++)
    {
        GRBLinExpr ptot = 0;
        for (int w = 0; w < 1; w++)
        {
            ptot += Sets[w][p]*t_var[w];
        }
         
        ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que //oferece capacidades para volverse string
        cname << "Battery " << p;
        master.addConstr(ptot <= 1, cname.str());
    }
    Sets.clear();
     
    //save the constraints in the array named c
    master.update();
     
    return master;
}
 
vector<double> Aux_MultiRole(int Nodes, int Targets,int Nroles,  double dual[], Mat &A1,vector<double> &nrj,double *ObjAux){
    //Declaration of environments and variables for the model
    GRBEnv* env = 0;
    GRBVar* E = 0;
    GRBVar** y = 0;
    GRBVar** x = 0;
    vector<double> Column;
     
     
    int count=0;
     
    double FObj;
    try{
         
        // Model
        env = new GRBEnv(); //creation of model environment
        GRBModel model = GRBModel(*env); //Declaration ofmodel
        model.set(GRB_StringAttr_ModelName, "Aux_MultiRole"); //Model Name Declaration
        model.getEnv().set(GRB_DoubleParam_Heuristics ,0.05);
        //model.getEnv().set(GRB_DoubleParam_IntFeasTol,1e-9);
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
         
        //***********************************************************
         
        // Total energy consumption decision variable
         
         
        E = model.addVars(Nodes, GRB_CONTINUOUS);
        model.update();
        for (int i = 0; i < Nodes; i++)
        {
            ostringstream vname;
            vname << "E_" << i;
            if(i>0){
                E[i].set(GRB_DoubleAttr_Obj, dual[i-1]);
                E[i].set(GRB_StringAttr_VarName, vname.str());
                E[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
            }
             
            if(i==0){
                E[i].set(GRB_DoubleAttr_Obj, 0);
                E[i].set(GRB_StringAttr_VarName, vname.str());
                E[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
            }
             
        }
         
        //************************** ******* *************************
         
        //************************** ******* *************************
        //Creating flow variables
         
        x = new GRBVar* [Nodes+Targets];
        for (int w = 0; w < Nodes+Targets; ++w)
        {
             
            x[w] = model.addVars(Nodes+Targets,GRB_INTEGER);
            count++;
            model.update();
            for (int p = 0; p < Nodes+Targets; ++p)
            {
                ostringstream vname;
                vname << "X_" << w << "_" << p;
                x[w][p].set(GRB_DoubleAttr_Obj, 0);
                x[w][p].set(GRB_StringAttr_VarName, vname.str());
            }
        }
         
         
        //***********************************************************
        //Creating Activation mode variables
        y = new GRBVar* [Nodes];
        for (int w = 0; w < Nodes; ++w) //Three is for the number of activation modes
        {
             
            y[w] = model.addVars(nrj.size(),GRB_BINARY);
            count++;
            model.update();
            for (int p = 0; p < nrj.size(); ++p)
            {
                ostringstream vname;
                vname << "y_" << w << "-" << p;
                y[w][p].set(GRB_DoubleAttr_Obj, 0);
                y[w][p].set(GRB_StringAttr_VarName, vname.str());
            }
        }
         
        model.update(); //Updating model to include variables
         
        //*************************************************** Including Balance Constraints  *************************************************
         
        // Flow constraints concerning to nodes
         
        for (int v = 0; v < Nodes+Targets; v++)
        {
            GRBLinExpr Act_Mode = 0;
            if (v==0) {
                for (int u = 0; u < Nodes+Targets; u++)
                {
                    Act_Mode += x[v][u]*A1[v][u];
                    Act_Mode += x[u][v]*A1[u][v]*(-1);
                     
                }
                ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
                cname << "Balance_" << v;
                model.addConstr(Act_Mode == -Targets, cname.str());
            }
            else{if (v>0 and v<Nodes) {
                for (int u = 0; u < Nodes+Targets; u++)
                {
                    Act_Mode += x[v][u]*A1[v][u];
                    Act_Mode += x[u][v]*A1[u][v]*(-1);
                     
                }
                ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
                cname << "Balance_" << v;
                model.addConstr(Act_Mode == 0, cname.str());
            }
            else{
                 
                for (int u = 0; u < Nodes+Targets; u++)
                {
                    Act_Mode += x[v][u]*A1[v][u];
                    Act_Mode += x[u][v]*A1[u][v]*(-1);
                     
                }
                ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
                cname << "Balance" << v;
                model.addConstr(Act_Mode == 1, cname.str());
                 
            }
            }
             
        }
         
         
         
         
        //*************************************************** Including Energy Level Constraint  *************************************************
         
        // Role Allocation to Sensors
        for (int v = 1;v < Nodes; v++)
        {
            GRBLinExpr Role_Alloc = 0;
            for (int l = 0; l < 3; l++){
                Role_Alloc += y[v][l];
            }
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "Role_Alloc_" << v;
            model.addConstr(Role_Alloc == 1, cname.str());
        }
         
        //Computing energy consumption for the selected activation mode
        for (int v = 1;v < Nodes; v++)
        {
            GRBLinExpr Energy_Alloc = 0;
            for (int l = 0; l < 3; l++){
                 
                Energy_Alloc += y[v][l]*nrj[l];
            }
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "Energy_Alloc" << v;
            model.addConstr(Energy_Alloc == E[v], cname.str());
        }
         
         
        for (int v = Nodes;v < Nodes+Targets; v++)
        {
            GRBLinExpr User_Cuts = 0;
            for (int u = 1; u < Nodes; u++){
                User_Cuts += y[u][2]*A1[v][u];
            }
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "User_Cuts_" << v;
            model.addConstr(User_Cuts >= 1, cname.str());
        }
         
        for (int v = 0;v < 1; v++)
        {
            GRBLinExpr User_Cuts = 0;
            for (int u = 1; u < Nodes; u++){
                User_Cuts += (y[u][2]+y[u][1])*A1[u][v];
            }
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "User_Cuts_" << v;
            model.addConstr(User_Cuts >= 1, cname.str());
        }
        //*************************************************** Coverage Constraints  *************************************************
        /*
         //Including coverage Constraints
         for (int u = Nodes;u < Nodes+Targets; u++)
         {
         GRBLinExpr Cov_Constraint = 0;
         for (int v = 1; v < Nodes; v++){
         Cov_Constraint += y[v][2]*A1[u][v];
         }
          
         ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
         cname << "Cov_Constraint" << u;
         model.addConstr(Cov_Constraint >= 1, cname.str());
         }*/
         
         
        //*************************************************** Flow Bounds based on connectivit and coverage  *****************************************
         
        // Role Allocation to Sensors
        for (int v = 0;v < Nodes; v++)
        {
            for (int u = 0; u < Nodes; u++){
                GRBLinExpr FlowBounds = 0;
                for (int l=1; l<3; l++) {
                    FlowBounds += y[v][l];
                }
                ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
                cname << "Bound_Sensor_sensor" << u<<"_"<<v;
                model.addConstr(A1[u][v]*x[u][v] <= FlowBounds*Targets, cname.str());
            }
        }
         
        for (int v = 0;v < Nodes; v++)
        {
            for (int u = 0; u < Nodes; u++){
                GRBLinExpr FlowBounds = 0;
                for (int l=1; l<3; l++) {
                    FlowBounds += y[v][l];
                }
                ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
                cname << "Bound_Sensor_sensor" << u<<"_"<<v;
                model.addConstr(A1[v][u]*x[v][u] <= FlowBounds*Targets, cname.str());
            }
        }
         
        for (int v = 0;v < Nodes; v++)
        {
            for (int u = Nodes; u < Nodes+Targets; u++){
                GRBLinExpr FlowBounds = 0;
                FlowBounds += y[v][2];
                ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
                cname << "Bound_Target_sensor" << u<<"_"<<v;
                model.addConstr(A1[u][v]*x[u][v] <= FlowBounds, cname.str());
            }
        }
         
         
        model.update();
         
        //model.write("/Users/facastanog/Desktop/LP_Auxiliar.lp");
         
         
        //Optimizing model to return column
        model.optimize();
         
         
        //Returning Objective fucntioj
        FObj=model.get(GRB_DoubleAttr_ObjVal);
         
         
        //Returning Column containing energy consumption information
        for(int i=1;i<Nodes;i++){
            Column.push_back(E[i].get(GRB_DoubleAttr_X));
            //~ cout<<E[i].get(GRB_DoubleAttr_X)<<" ";
        }
        //cout<<endl;
         
        /*for(int i=0;i<Nodes+Targets;i++){
         for(int j=0;j<Nodes+Targets;j++){
         if (x[i][j].get(GRB_DoubleAttr_X)>0) {
         cout<<i<<" "<<j<<" "<<x[i][j].get(GRB_DoubleAttr_X)<<endl;
         }
          
         }
         }
         */
        *ObjAux=FObj;
         
         
    }
    catch(...){
        cout<<"Error in Auxiliary problem creation";
         
    }
    delete[] y;
    delete[] x;
    return Column;
     
}
 
Mat Con_Comp(Mat &Adj, Vec &Solution){
    //Function used to store the connected components present in the solution by using an integer representation that could be used by floyd Warshall to identify the shortest path to join those components
     
    Vec Labels,Comp;
    Mat ConnComp;
    for (int i=0; i<Solution.size(); i++) {
        Labels.push_back(0);
    }
     
    Comp.push_back(0);
    Labels[0]=1;
    DFS(Adj, Solution, Labels, 0, Comp);
    ConnComp.push_back(Comp);
    Comp.clear();
     
    for(int i=0; i<Solution.size(); i++) {
        //A DFS is performed recording in "Labels" the nodes already visited by using pointers
        if (Solution[i]>0.0000001 and Labels[i]==0) {
            DFS(Adj, Solution, Labels, i, Comp);
            ConnComp.push_back(Comp);
            Comp.clear();
        }
    }
    return ConnComp;
     
}
 
 
 
void DFS(Mat &Adj, Vec const &Solution, Vec &Visited, int cur_node, Vec &Comp){
    int cont=0;
    while (cont<Solution.size()){
        if ((Adj[cur_node][cont]==1 or Adj[cont][cur_node]==1) and Solution[cont]>0.0000001 and Visited[cont]==0){
            Visited[cont]=1;
            Comp.push_back(cont);
            DFS(Adj, Solution, Visited, cont, Comp);
        }
        cont+=1;
         
    }
}

int Check_cov(Mat &Adj, Vec const &Solution, int T, int S){
    int cont=0;
    int Tot_Cov=0;
    Vec cover;
    for (int i=0; i<T; i++) {
        cover.push_back(0);
    }
    while (cont<Solution.size()){
        for (int i=0; i<T; i++) {
            if(cover[i]==0 and Adj[S+1+i][Solution[cont]]==1){
                cover[i]=1;
                Tot_Cov+=1;
                 
            }
        }
        cont+=1;
    }
    return Tot_Cov;
}
