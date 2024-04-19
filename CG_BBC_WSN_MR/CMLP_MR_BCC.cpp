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
//g++ -Ofast -m64 -g -o main_CG_BBC.out CG_BBC_MR.cpp -I/opt/gurobi562/linux64/include -L/opt/gurobi562/linux64/lib -lgurobi_c++ -L/opt/gurobi562/linux64/lib -lgurobi56 -lpthread -lm

 
//Data structure definitions
typedef vector<int> Vec;
typedef vector<Vec> Mat;
typedef vector<Mat> Array3D;
 
//initialization
Vec Random_Sol(Mat &Adj,int Nodes, int T, int *Act_Sen, double alpha);
void RND_DFS(Mat &Adj, Vec &Sol,Vec &Cob, int Nodes, int T, int Cur_Node,int *CovTotal, double alpha ); //Generate connected random initial solution
 
 
//Declarations
double metodos(string dir,   int method,double Rs,double Rc);
double Column_Generation_ILP(string dir,  int method,double Rs,double Rc);
double Column_Generation_CP(string dir,  int method,double Rs,double Rc);
double Stab_CG_CP(string dir,  int method,double Rs,double Rc);
double Neame_Generation_CP(string dir,  int method,double Rs,double Rc);
 
int ConnBin(double x1, double x2, double y1, double y2,double Rs);
vector<double> Aux_MultiRole(int Nodes, int Targets,int Nroles,  double dual[], Mat &A1,vector<double> &nrj,double *ObjAux); //Function returning the best column found
 
 
//Benders declarations
vector<double> Benders_Flow(Mat &A1, Mat &y, int Nodes, int Targets, double *Bound,int *status);
vector<double> Benders_Algorithm(int Nodes, int Targets,int Nroles,  double dual[], Mat &A1,vector<double> &nrj,double *ObjAux);
vector<double> Benders_Algorithm_Complete(int Nodes, int Targets,int Nroles,  double dual[], Mat &A1,vector<double> &nrj,double *ObjAux);
 
//Functions used to create the models
GRBModel Gen_Master(int S,vector<double> &nrj);
GRBModel Gen_Aux_MultiRole(int Nodes, int Targets,int Nroles, Mat &A1,vector<double> &nrj);
GRBModel Gen_Cutting_planes(int Nodes, int Targets,int Nroles,   Mat const &A1,vector<double> &nrj);
 
//Functions used to identify connected components
Mat Con_Comp(Mat &Adj, Vec &Solution);
void DFS(Mat &Adj, Vec const &Solution, Vec &labels, int cur_node, Vec &Comp);
 
//Callback class definition
vector<double> Cutting_planes(int Nodes, int Targets,int Nroles,  double dual[], Mat &A1,vector<double> &nrj,double *ObjAux);
 
 
//Callback related stuff
vector< vector<int> > NewCols;
int Check_cov(Mat &Adj, Vec const &Solution, int T, int S);
 
//global model used to manage BSP
GRBEnv* envBenders = new GRBEnv(); //creation of model environment
GRBModel BSP = GRBModel(*envBenders);
 
 
class connectivity: public GRBCallback
{
public:
    GRBVar** vars;
    int n_nodes;
    int nrj_lev;
    int targets;
    vector< vector<int> > Adj;
    vector< vector<int> > Cols;
    connectivity(GRBVar** Y_vars, int xn, int xl,int T,vector< vector<int> > const &A1) {
        vars = Y_vars;
        n_nodes = xn;
        nrj_lev=xl;
        Adj=A1;
        targets=T;
    }
protected:
    void callback() {
        try {
            vector<int> x;
            vector<int> members;
            Mat CC;
            int coverage=0;
            switch (where) {
               case GRB_CB_MIPSOL:
                    // Found an integer feasible solution - does it visit every node?
                     
                    for (int i = 0; i < n_nodes; i++){
                        //The numbers one and two are used to denote the index of the Y solution according to the representation in the solutions
                        x.push_back(round(getSolution(vars[i][1])+2*getSolution(vars[i][2])));
                    }
                    //x[0]=1;
                    CC=Con_Comp(Adj, x);
                    for (int i=0; i<CC[0].size(); i++) {
                        if (x[CC[0][i]]==2) {
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
                        //Add connectiivity elimination constraint
                        for (int i=1; i<CC.size(); i++) {
                            GRBLinExpr left = 0;
                            GRBLinExpr right = 0;
                            for (int j=0; j<CC[i].size(); j++) {
                                left+=vars[CC[i][j]][1]+vars[CC[i][j]][2];
                                for (int k=0; k<n_nodes; k++) {
                                    if ((Adj[CC[i][j]][k]>0 or Adj[k][CC[i][j]]>0) and (x[k]==0)) {
                                        right+=vars[k][1]+vars[k][2];
                                    }
                                }
                            }
                            addLazy(left<=CC[i].size()*right);
                            //cout<<left<<"<="<<CC[i].size()<<"*"<<right<<endl;
                        }
                         
                    }
                    break;
  
 
            }
             
        } catch (GRBException e) {
            cout << "Error number: " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) {
            cout << "Error during callback" << endl;
        }
    }
};
 
class connectivity_Aggressive: public GRBCallback
{
public:
    GRBVar** vars;
    int n_nodes;
    int nrj_lev;
    int targets;
    vector< vector<int> > Adj;
    vector< vector<int> > Cols;
    connectivity_Aggressive(GRBVar** Y_vars, int xn, int xl,int T,vector< vector<int> > const &A1) {
        vars = Y_vars;
        n_nodes = xn;
        nrj_lev=xl;
        Adj=A1;
        targets=T;
    }
protected:
    void callback() {
        try {
            vector<int> x;
            vector<int> members;
            Mat CC;
            int coverage=0;
            switch (where) {
                case GRB_CB_MIPSOL:
                    // Found an integer feasible solution - does it visit every node?
                     
                    for (int i = 0; i < n_nodes; i++){
                        //The numbers one and two are used to denote the index of the Y solution according to the representation in the solutions
                        x.push_back(round(getSolution(vars[i][1])+2*getSolution(vars[i][2])));
                    }
                    //x[0]=1;
                    CC=Con_Comp(Adj, x);
                    for (int i=0; i<CC[0].size(); i++) {
                        if (x[CC[0][i]]==2) {
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
                        //Add connectiivity elimination constraint
                        for (int i=1; i<CC.size(); i++) {
                            GRBLinExpr left = 0;
                            GRBLinExpr right = 0;
                            for (int j=0; j<CC[i].size(); j++) {
                                left+=vars[CC[i][j]][1]+vars[CC[i][j]][2];
                                for (int k=0; k<n_nodes; k++) {
                                    if ((Adj[CC[i][j]][k]>0 or Adj[k][CC[i][j]]>0) and (x[k]==0)) {
                                        right+=vars[k][1]+vars[k][2];
                                    }
                                }
                            }
                            addLazy(left<=CC[i].size()*right);
                            //cout<<left<<"<="<<CC[i].size()<<"*"<<right<<endl;
                        }
                         
                    }
                    break;
                case GRB_CB_MIPNODE:
                    //cout<<"here I am "<<getIntInfo(GRB_CB_MIPNODE_STATUS)<<" "<<GRB_OPTIMAL<<endl;
                    if (getIntInfo(GRB_CB_MIPNODE_STATUS)==GRB_OPTIMAL) {
						x.clear();
                        // Found an integer feasible solution - does it visit every node?
                        for (int i = 0; i < n_nodes; i++){
                            //The numbers one and two are used to denote the index of the Y solution according to the representation in the solutions
                            x.push_back(round(getNodeRel(vars[i][1]))+2*round(getNodeRel(vars[i][2])));
                        }
                        x[0]=1;
                        CC=Con_Comp(Adj, x);
						for (int i=0; i<CC[0].size(); i++) {
							if (x[CC[0][i]]==2) {
								members.push_back(CC[0][i]);
							}
						}
						coverage=Check_cov(Adj, members,targets, n_nodes-1);
						if (CC.size()>1 and coverage<targets) {
                            // Add connectiivity elimination constraint
								for (int i=1; i<CC.size(); i++) {
									GRBLinExpr left = 0;
									GRBLinExpr right =0;
									for (int j=0; j<CC[i].size(); j++) {
										left+=vars[CC[i][j]][1]+vars[CC[i][j]][2];
										for (int k=0; k<n_nodes; k++) {
											if ((Adj[CC[i][j]][k]>0 or Adj[k][CC[i][j]]>0) and (x[k]==0)) {
												right+=vars[k][1]+vars[k][2];
											}
										}
									}
									addLazy(left<=CC[i].size()*right);
									//cout<<left<<"<="<<CC[i].size()<<"*"<<right<<endl;
								}
								 
						}
						
                    }
                    break;
                     
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
		for (int size=1; size<=1; size++) {
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
						try{
							MLP=Column_Generation_CP(dir,3,RadS[m],RadC[m]);
						}
							catch (GRBException e) {
							cout << "weirdest: " << e.getErrorCode() << endl;
							cout << e.getMessage() << endl;
						} catch (...) {
							cout << "This is weird dude" << endl;
						}
					}
				}
			}
     }
     }//*/
     
    time_t start,end;
    time (&start);
    dir="/home/fabian/Documents/MLP/Instance_1_200_15_3.txt";
    MLP=Column_Generation_CP(dir,2,125,125);
    //MLP=Stab_CG_CP(dir,2,125,125);
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
     
     
    vector<double> nrj;
    nrj.push_back(0);
    nrj.push_back(0.8);
    nrj.push_back(1.0);
     
     
    //Preparing the creation of the Model
    int ncols=0; //Variable to store the number of generated columns
    GRBModel MP=Gen_Master(S,nrj);
     
    //=================================================================
    GRBEnv* envPS = 0;
    envPS = new GRBEnv(); //creation of model environment
    GRBVar** Y = 0; //variables concerning to the use of the sensors and the related state
    Y= new GRBVar* [S+1];
     
    GRBModel PS = GRBModel(*envPS); //Declaration ofmodel
    PS.set(GRB_StringAttr_ModelName, "Benders_Master"); //Model Name Declaration
    PS.getEnv().set(GRB_IntParam_OutputFlag, 0);
    PS.getEnv().set(GRB_IntParam_LazyConstraints, 1);
    PS.getEnv().set(GRB_DoubleParam_Heuristics,0.0);
    PS.getEnv().set(GRB_IntParam_MIPFocus,1);
     
    // The objective is to maximize the total fixed and variable costs
    PS.set(GRB_IntAttr_ModelSense, 1);
    PS.getEnv().set(GRB_DoubleParam_TimeLimit, 10);
    //model.getEnv().set(GRB_IntParam_Presolve,0);
     
     
    for (int w = 0; w < S+1; ++w) //Three is for the number of activation modes
    {
        Y[w] = PS.addVars(nrj.size(),GRB_BINARY);
        PS.update();
        for (int p = 0; p < nrj.size(); ++p)
        {
            ostringstream vname;
            vname << "Y_" << w << "-" << p;
            Y[w][p].set(GRB_DoubleAttr_Obj, 0);
            Y[w][p].set(GRB_StringAttr_VarName, vname.str());
        }
    }
     
    PS.update();
    //Adding initial cuts
    // Role Allocation to Sensors
    for (int v = 0;v < S+1; v++)
    {
        GRBLinExpr Role_Alloc = 0;
        for (int l = 0; l < 3; l++){
            Role_Alloc += Y[v][l];
        }
        ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
        cname << "Role_Alloc_" << v;
        PS.addConstr(Role_Alloc == 1, cname.str());
    }
    for (int v = 0;v < 1; v++)
    {
        GRBLinExpr User_Cuts = 0;
        for (int u = 1; u < S+1; u++){
            User_Cuts += (Y[u][2]+Y[u][1])*A1[u][v];
        }
        ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
        cname << "User_Cuts_" << v;
        PS.addConstr(User_Cuts >= 1, cname.str());
    }


    for (int v = S+1;v < S+1+T; v++)
    {
        GRBLinExpr User_Cuts = 0;
        for (int u = 1; u < S+1; u++){
            User_Cuts += Y[u][2]*A1[v][u];
        }
        ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
        cname << "User_Cuts_" << v;
        PS.addConstr(User_Cuts >= 1, cname.str());
    }
     // /*    
    for (int v = 1;v < S+1; v++)
    {
        if(A1[v][0]<=0){
            GRBLinExpr User_Cuts = 0;
            for (int u = 1; u < S+1; u++){
                User_Cuts += (Y[u][2]+Y[u][1])*A1[u][v];
            }
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "User_Cuts_" << v;
            PS.addConstr(User_Cuts >= Y[v][2]+Y[v][1], cname.str());
        }
         
    }//*/
       
    GRBLinExpr User_Cuts = 0;
    for (int v = 1;v < S+1; v++)
    {
        User_Cuts += (Y[v][1]+Y[v][2]);
    }
    ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
    cname << "Max_relays";
    int RHS=30,Remo_ctr=0;
    PS.addConstr(User_Cuts <= RHS, cname.str());
   
    PS.update();
    //The model is solved calling the "callback" used to add the connectivity constraints
    //model.write("/Users/Fabian/Dropbox/Imprimid.lp");
     
    //solution pool: is used to sttore the solutions found along the CG process
    connectivity cb=connectivity(Y, S+1,nrj.size(), T,A1);
    connectivity_Aggressive cb2=connectivity_Aggressive(Y, S+1,nrj.size(), T,A1);
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

	PS.setCallback(&cb);
	PS.getEnv().set(GRB_DoubleParam_Cutoff,1);
    while (true) {
        iters+=1;
        MP.optimize();
        time(&end);
        cout<<" Lifetime: "<<MP.get(GRB_DoubleAttr_ObjVal)<<" Time: "<<difftime(end, start)<<" n cols: " << ncols<<" n iters: " << iters<<endl;
         
        //This is a subroutine used to modify the objective function of the PS according to the duals found throught the MP!
        for(int i=0;i<S;i++){
            for (int j=1; j<nrj.size(); j++) {
                ostringstream sname;
                //The constraint i is related to the variables of the Nodes i+1
                sname << "Y_" << i+1 <<"-"<<j;
                Y[i+1][j].set(GRB_DoubleAttr_Obj,c[i].get(GRB_DoubleAttr_Pi)*nrj[j]);
            }
        }
        //PS.update();
        //PS.write("/Users/Fabian/Desktop/Salida.lp");
        NewCols.clear();
        PS.optimize();
        time(&end);
        if(PS.get(GRB_IntAttr_Status)!=3 and 1-PS.get(GRB_DoubleAttr_ObjVal)>0.000001 and  (PS.get(GRB_IntAttr_Status)==2 or PS.get(GRB_IntAttr_Status)==13) and difftime(end, start)<3600){
            
          for (int p=NewCols.size()-1; p<NewCols.size(); p++) {
                ncols++;
                //cout <<"Adding Col "<<ncols<<" Lifetime:"<<MP.get(GRB_DoubleAttr_ObjVal)<<endl;
                GRBColumn Col;
                //The columns, in the RMP, are created with the energy of the sensors exactly as they are stacked in the E_x vector, whic starts with the first sensor
                double E_tot=0;
                for (int i=0; i<S; i++) {
                    //The constraint i is related to the variables of the Nodes i+1
                    E_tot=nrj[NewCols[p][i+1]];
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
            if ((PS.get(GRB_IntAttr_Status)!=GRB_OPTIMAL or PS.get(GRB_IntAttr_Status)==3) and  difftime(end, start)<3600) {
				cont_sense+=1;
				PS.getEnv().set(GRB_DoubleParam_TimeLimit, 10*(cont_sense+1));
				
				if(NewCols.size()>0){
					cont_sense=0;
					for (int p=0; p<NewCols.size(); p++) {
                    ncols++;
                    //cout <<"Adding Col "<<ncols<<" Lifetime:"<<MP.get(GRB_DoubleAttr_ObjVal)<<endl;
                    GRBColumn Col;
                    //The columns, in the RMP, are created with the energy of the sensors exactly as they are stacked in the E_x vector, whic starts with the first sensor
                    double E_tot=0;
                    for (int i=0; i<S; i++) {
                        //The constraint i is related to the variables of the Nodes i+1
                        E_tot=nrj[NewCols[p][i+1]];
                        //cout<<E_tot<<" ";
                        if (E_tot>0)
                            Col.addTerm(E_tot, c[i]);
                    }                     
                    //cout<<endl;
                    MP.addVar(0.0, GRB_INFINITY, 1, GRB_CONTINUOUS,Col,"t"+ncols);
                }
                if (Remo_ctr==0) {
                    cout<<"Aggressive"<<endl;
                    PS.remove(Ps_ctr[PS.get(GRB_IntAttr_NumConstrs)-1]);
                    Remo_ctr=1;
                }
			}
                
            }
            else{
				FILE_Exact_Results.open("/home/fabian/Documents/CG_BBC_WSN_MR/Res_CG_BBC_MR_300.txt",fstream::app);
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
    nrj.push_back(0.8);
    nrj.push_back(1.0);
     
    //Preparing the creation of the Model
    int ncols=0; //Variable to store the number of generated columns
    GRBModel MP=Gen_Master(S,nrj);
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
                PS.getVarByName(sname.str()).set(GRB_DoubleAttr_Obj,c[i].get(GRB_DoubleAttr_Pi)*nrj[j]);
            }
        }
        //PS.update();
        //PS.write("/Users/Fabian/Desktop/Salida.lp");
        PS.optimize();
         
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
                    FILE_Exact_Results.open("//Users/Fabian/Desktop/Res_CG+ILP+Stab.txt",fstream::app);
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
        int RHS=15;
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
 
GRBModel Gen_Master(int S,vector<double> &nrj){
    //Environment and model related variables
    GRBEnv* envMP = 0;
    GRBVar* t_var = 0;
     
    //Initialization
    vector< vector<double> > Sets;
    vector<double> ColIni;
    for (int i=0; i<S; i++) {
        ColIni.push_back(nrj[2]);
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
 
vector<double> Cutting_planes(int Nodes, int Targets,int Nroles,  double dual[], Mat &A1,vector<double> &nrj,double *ObjAux){
     
    GRBEnv* env = 0;
    env = new GRBEnv(); //creation of model environment
    GRBVar** Y = 0; //variables concerning to the use of the sensors and the related state
    Y= new GRBVar* [Nodes];
     
    GRBModel model = GRBModel(*env); //Declaration ofmodel
    model.set(GRB_StringAttr_ModelName, "Benders_Master"); //Model Name Declaration
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
    //model.getEnv().set(GRB_DoubleParam_Heuristics,0);
    // The objective is to maximize the total fixed and variable costs
    model.set(GRB_IntAttr_ModelSense, 1);
    //model.getEnv().set(GRB_IntParam_Presolve,0);
     
     
    for (int w = 0; w < Nodes; ++w) //Three is for the number of activation modes
    {
        if(w==0){
             
             
            Y[w] = model.addVars(nrj.size(),GRB_BINARY);
            model.update();
            for (int p = 0; p < nrj.size(); ++p)
            {
                ostringstream vname;
                vname << "Y_" << w << "-" << p;
                Y[w][p].set(GRB_DoubleAttr_Obj, 0);
                Y[w][p].set(GRB_StringAttr_VarName, vname.str());
            }
        }
        if(w>0){
            Y[w] = model.addVars(nrj.size(),GRB_BINARY);
            model.update();
            for (int p = 0; p < nrj.size(); ++p)
            {
                ostringstream vname;
                vname << "Y_" << w << "-" << p;
                Y[w][p].set(GRB_DoubleAttr_Obj, nrj[p]*dual[w-1]);
                Y[w][p].set(GRB_StringAttr_VarName, vname.str());
            }
        }
    }
     
    model.update();
    //Adding initial cuts
    // Role Allocation to Sensors
    for (int v = 0;v < Nodes; v++)
    {
        GRBLinExpr Role_Alloc = 0;
        for (int l = 0; l < 3; l++){
            Role_Alloc += Y[v][l];
        }
        ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
        cname << "Role_Alloc_" << v;
        model.addConstr(Role_Alloc == 1, cname.str());
    }
    for (int v = Nodes;v < Nodes+Targets; v++)
    {
        GRBLinExpr User_Cuts = 0;
        for (int u = 1; u < Nodes; u++){
            User_Cuts += Y[u][2]*A1[v][u];
        }
        ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
        cname << "User_Cuts_" << v;
        model.addConstr(User_Cuts >= 1, cname.str());
    }
     
    for (int v = 0;v < 1; v++)
    {
        GRBLinExpr User_Cuts = 0;
        for (int u = 1; u < Nodes; u++){
            User_Cuts += (Y[u][2]+Y[u][1])*A1[u][v];
        }
        ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
        cname << "User_Cuts_" << v;
        model.addConstr(User_Cuts >= 1, cname.str());
    }
     
    for (int v = 1;v < Nodes; v++)
    {
        if(A1[v][0]<=0){
            GRBLinExpr User_Cuts = 0;
            for (int u = 1; u < Nodes; u++){
                User_Cuts += (Y[u][2]+Y[u][1])*A1[u][v];
            }
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "User_Cuts_" << v;
            model.addConstr(User_Cuts >= Y[v][2]+Y[v][1], cname.str());
        }
         
    }
     
    model.update();
    //The model is solved calling the "callback" used to add the connectivity constraints
    //model.write("/Users/Fabian/Dropbox/Imprimid.lp");
    //solution pool: is used to sttore the solutions found along the CG process
    vector< vector<int> > SolPool;
    connectivity cb=connectivity(Y, Nodes, nrj.size(),Targets, A1);
    model.setCallback(&cb);
    model.optimize();
     
     
     
    //Returning Column containing energy consumption information
    double E;
    vector<double> Column;
    for(int i=1;i<Nodes;i++){
        E=0;
        for (int l=0; l<nrj.size(); l++) {
            E+=Y[i][l].get(GRB_DoubleAttr_X)*nrj[l];
        }
        Column.push_back(E);
        // cout<<E<<", ";
    }
    //cout<<endl;
    *ObjAux=model.get(GRB_DoubleAttr_ObjVal);
    //cout<<"Obj_Sol: "<<model.get(GRB_DoubleAttr_ObjVal)<<endl;
    for (int i=0; i<Nodes; i++) {
        delete[] Y[i];
    }
    delete[] Y;
    delete env;
    return Column;
     
}
 
vector<double> Benders_Algorithm_Complete(int Nodes, int Targets,int Nroles,  double dual[], Mat &A1,vector<double> &nrj,double *ObjAux){
    //Creating master model
    double Lower_Bound=-100000000, Upper_Bound=0;
     
    GRBEnv* env = 0;
    env = new GRBEnv(); //creation of model environment
    GRBVar** Y = 0; //variables concerning to the use of the sensors and the related state
    Y= new GRBVar* [Nodes];
     
    GRBModel model = GRBModel(*env); //Declaration ofmodel
    model.set(GRB_StringAttr_ModelName, "Benders_Master"); //Model Name Declaration
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    //model.getEnv().set(GRB_DoubleParam_Heuristics,0);
    // The objective is to maximize the total fixed and variable costs
    model.set(GRB_IntAttr_ModelSense, 1);
    //model.getEnv().set(GRB_IntParam_Presolve,0);
     
     
    Mat y;
    Vec Fila;
    int status=0;
    vector<double> ray;
    int tracker=0;
    double New_Upper_Bound,Bound;
    int Ncuts=0;
     
     
    for (int w = 0; w < Nodes; ++w) //Three is for the number of activation modes
    {
        if(w==0){
             
             
            Y[w] = model.addVars(nrj.size(),GRB_BINARY);
            model.update();
            for (int p = 0; p < nrj.size(); ++p)
            {
                ostringstream vname;
                vname << "Y_" << w << "-" << p;
                Y[w][p].set(GRB_DoubleAttr_Obj, 0);
                Y[w][p].set(GRB_StringAttr_VarName, vname.str());
            }
        }
        if(w>0){
            Y[w] = model.addVars(nrj.size(),GRB_BINARY);
            model.update();
            for (int p = 0; p < nrj.size(); ++p)
            {
                ostringstream vname;
                vname << "Y_" << w << "-" << p;
                Y[w][p].set(GRB_DoubleAttr_Obj, nrj[p]*dual[w-1]);
                Y[w][p].set(GRB_StringAttr_VarName, vname.str());
            }
        }
    }
     
    model.update();
     
    //Adding initial cuts
     
     
    // Role Allocation to Sensors
    for (int v = 0;v < Nodes; v++)
    {
        GRBLinExpr Role_Alloc = 0;
        for (int l = 0; l < 3; l++){
            Role_Alloc += Y[v][l];
        }
        ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
        cname << "Role_Alloc_" << v;
        model.addConstr(Role_Alloc == 1, cname.str());
    }
    for (int v = Nodes;v < Nodes+Targets; v++)
    {
        GRBLinExpr User_Cuts = 0;
        for (int u = 1; u < Nodes; u++){
            User_Cuts += Y[u][2]*A1[v][u];
        }
        ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
        cname << "User_Cuts_" << v;
        model.addConstr(User_Cuts >= 1, cname.str());
    }
     
    for (int v = 0;v < 1; v++)
    {
        GRBLinExpr User_Cuts = 0;
        for (int u = 1; u < Nodes; u++){
            User_Cuts += (Y[u][2]+Y[u][1])*A1[u][v];
        }
        ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
        cname << "User_Cuts_" << v;
        model.addConstr(User_Cuts >= 1, cname.str());
    }
     
     
     
    //Routines used to store initial solution and to compute initial upper bounds
     
    for(int i=0;i<Nodes;i++){
        for (int j=0; j<nrj.size(); j++) {
            if (j<nrj.size()-1) {
                Fila.push_back(0);
            }
            else{
                Fila.push_back(1);
            }
        }
        y.push_back(Fila);
        Fila.clear();
    }
     
    for(int i=1;i<Nodes;i++){
        for (int l=0; l<nrj.size(); l++) {
            Upper_Bound+=nrj[l]*dual[i-1]*y[i][l];
        }
    }
     
    // =============================================================Creating auxiliary model================================================================
     
    //cout<<"creating Aux_model"<<endl;
    //Declarations
    GRBEnv* env_aux = 0;
    GRBVar** x = 0;
    GRBConstr* c = 0;
    x = new GRBVar* [Nodes+Targets];
    env_aux = new GRBEnv(); //creation of model environment
    GRBModel BSP = GRBModel(*env_aux); //Declaration ofmodel
    BSP.set(GRB_IntAttr_ModelSense, 1);
     
    //Setting GUROBI parameters
    BSP.set(GRB_StringAttr_ModelName, "Benders_Master"); //Model Name Declaration
    BSP.getEnv().set(GRB_IntParam_OutputFlag, 0); //Screen log
    BSP.getEnv().set(GRB_IntParam_Presolve,0);
    BSP.getEnv().set(GRB_DoubleParam_IntFeasTol,1e-9); //Feasibility tolerance
    BSP.getEnv().set(GRB_IntParam_InfUnbdInfo, 1); //Allows to obtain farkas dual information
     
    //************************** ******* *************************
    //Creating flow variables
     
    for (int w = 0; w < Nodes+Targets; ++w)
    {
        x[w] = BSP.addVars(Nodes+Targets,GRB_CONTINUOUS);
        BSP.update();
        for (int p = 0; p < Nodes+Targets; ++p)
        {
            ostringstream vname;
            vname << "X_" << w << "_" << p;
            x[w][p].set(GRB_DoubleAttr_Obj, 0);
            x[w][p].set(GRB_StringAttr_VarName, vname.str());
        }
    }
     
    BSP.update(); //Updating model to include variables
     
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
            BSP.addConstr(Act_Mode == -Targets, cname.str());
        }
        else{if (v>0 and v<Nodes) {
            for (int u = 0; u < Nodes+Targets; u++)
            {
                Act_Mode += x[v][u]*A1[v][u];
                Act_Mode += x[u][v]*A1[u][v]*(-1);
                 
            }
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "Balance_" << v;
            BSP.addConstr(Act_Mode == 0, cname.str());
        }
        else{
             
            for (int u = 0; u < Nodes+Targets; u++)
            {
                Act_Mode += x[v][u]*A1[v][u];
                Act_Mode += x[u][v]*A1[u][v]*(-1);
                 
            }
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "Balance" << v;
            BSP.addConstr(Act_Mode == 1, cname.str());
             
        }
        }
         
    }
     
     
     
     
    //*************************************************** Including Energy Level Constraint  *************************************************
     
    // Role Allocation to Sensors
    for (int v = 0;v < Nodes; v++)
    {
        for (int u = 0; u < Nodes; u++){
            GRBLinExpr FlowBounds = 0;
            for (int l=1; l<3; l++) {
                FlowBounds += y[v][l];
            }
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "Bound_Input_" << u<<"_"<<v;
            BSP.addConstr(A1[u][v]*x[u][v] <= FlowBounds*Targets, cname.str());
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
            cname << "Bound_Output_" << u<<"_"<<v;
            BSP.addConstr(A1[v][u]*x[v][u] <= FlowBounds*Targets, cname.str());
        }
    }
     
    for (int v = 0;v < Nodes; v++)
    {
        for (int u = Nodes; u < Nodes+Targets; u++){
            GRBLinExpr FlowBounds = 0;
            FlowBounds += y[v][2];
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "Bound_Target_Output_" << u<<"_"<<v;
            BSP.addConstr(A1[u][v]*x[u][v] <= FlowBounds, cname.str());
        }
    }
     
    BSP.update();
    //Auxiliary.write("/Users/facastanog/Desktop/Aux_Benders.lp");
    //cout<<"Finished creation"<<endl;
     
     
    // =================================================== Finishing creation of auxiliary model================================================================
     
    New_Upper_Bound=Upper_Bound;
    c=BSP.getConstrs();
    int iter=0;
     
     
    while(Upper_Bound>Lower_Bound){
         
        //cout<<"Lower_Bound iter "<<iter<<" : "<<Lower_Bound<<endl;
        iter+=1;
         
        //================================================Instructions performed to avoid AP to be created again   =============================================
        ray.clear();
         
        BSP.optimize();
        if (BSP.get(GRB_IntAttr_Status)==3) {
            status=3;
            New_Upper_Bound = BSP.get(GRB_DoubleAttr_ObjVal)+Lower_Bound;
             
             
            for(int i=0;i<BSP.get(GRB_IntAttr_NumConstrs);i++){
                ray.push_back(c[i].get(GRB_DoubleAttr_FarkasDual));
                 
            }
             
        }
        else{
            status=2;
        }
         
        //======================================================== Master problem related routines =============================================================
        if(status==2){
            Upper_Bound=New_Upper_Bound;
        }
        else{
             
            tracker=0;
            GRBLinExpr Feasibility_Cuts = 0;
            // Flow constraints concerning to nodes
             
            for (int v = 0; v < Nodes+Targets; v++)
            {
                if (v==0) {
                    Feasibility_Cuts += -Targets*ray[tracker];
                    tracker+=1;
                }
                else{if (v>0 and v<Nodes) {
                    Feasibility_Cuts += 0*ray[tracker];
                    tracker+=1;
                }
                else{
                    Feasibility_Cuts += ray[tracker];
                    tracker+=1;
                     
                }
                }
                 
            }
            // Role Allocation to Sensors
            for (int v = 0;v < Nodes; v++)
            {
                for (int u = 0; u < Nodes; u++){
                    GRBLinExpr FlowBounds = 0;
                    for (int l=1; l<3; l++) {
                        FlowBounds += Y[v][l];
                    }
                    Feasibility_Cuts+=ray[tracker]*FlowBounds*Targets;
                    tracker+=1;
                }
            }
            for (int v = 0;v < Nodes; v++)
            {
                for (int u = 0; u < Nodes; u++){
                    GRBLinExpr FlowBounds = 0;
                    for (int l=1; l<3; l++) {
                        FlowBounds += Y[v][l];
                    }
                    Feasibility_Cuts+=FlowBounds*Targets*ray[tracker];
                    tracker+=1;
                }
            }
            for (int v = 0;v < Nodes; v++)
            {
                for (int u = Nodes; u < Nodes+Targets; u++){
                    GRBLinExpr FlowBounds = 0;
                    FlowBounds += Y[v][2];
                    Feasibility_Cuts+= FlowBounds*ray[tracker];
                    tracker+=1;
                }
            }
             
            ostringstream cname;//objeto de tipo cname que es usado para guardar el nombre de las variables y que
            cname << "Feas_Cut_" << Ncuts;
            model.addConstr(Feasibility_Cuts >= 0, cname.str());
            model.update();
            //model.write("/Users/facastanog/Desktop/Master_Benders.lp");
        }
         
        model.optimize();
        //cout<<model.get(GRB_IntAttr_Status)<<endl;
        Lower_Bound=model.get(GRB_DoubleAttr_ObjVal);
         
        New_Upper_Bound=Lower_Bound;
         
        for (int i=0; i<Nodes; i++) {
            for (int l=0; l<nrj.size();l++){
                y[i][l]=Y[i][l].get(GRB_DoubleAttr_X);
            }
            //cout<<endl;
        }
        //================================================ Modifying RHS =============================================
        // Flow constraints concerning to nodes
        int Num_CTR=0;
        for (int v = 0; v < Nodes+Targets; v++)
        {
            GRBLinExpr Act_Mode = 0;
            if (v==0) {
                c[Num_CTR].set(GRB_DoubleAttr_RHS, -Targets);
                Num_CTR+=1;
            }
            else{if (v>0 and v<Nodes) {
                c[Num_CTR].set(GRB_DoubleAttr_RHS, 0);
                Num_CTR+=1;
            }
            else{
                c[Num_CTR].set(GRB_DoubleAttr_RHS, 1);
                Num_CTR+=1;
                 
            }
            }
             
        }
        //*************************************************** Including Energy Level Constraint  *************************************************
        // Role Allocation to Sensors
        for (int v = 0;v < Nodes; v++)
        {
            for (int u = 0; u < Nodes; u++){
                double Flows = 0;
                for (int l=1; l<3; l++) {
                    Flows += y[v][l];
                }
                c[Num_CTR].set(GRB_DoubleAttr_RHS, Flows*Targets);
                Num_CTR+=1;
            }
        }
        for (int v = 0;v < Nodes; v++)
        {
            for (int u = 0; u < Nodes; u++){
                double Flows = 0;
                for (int l=1; l<3; l++) {
                    Flows += y[v][l];
                }
                c[Num_CTR].set(GRB_DoubleAttr_RHS, Flows*Targets);
                Num_CTR+=1;
            }
        }
        for (int v = 0;v < Nodes; v++)
        {
            for (int u = Nodes; u < Nodes+Targets; u++){
                double Flows = 0;
                Flows += y[v][2];
                c[Num_CTR].set(GRB_DoubleAttr_RHS, Flows*Targets);
                Num_CTR+=1;
            }
        }
        BSP.update();
        //Auxiliary.write("/Users/Fabian/Desktop/Aux_Benders.lp");
        //cout<<"Finished creation"<<endl;
        //================================================ Finished RHS  =============================================
        Ncuts+=1;
    }//End_While
     
     
    //Returning Column containing energy consumption information
    double E;
    vector<double> Column;
    for(int i=1;i<Nodes;i++){
        E=0;
        for (int l=0; l<nrj.size(); l++) {
            E+=y[i][l]*nrj[l];
        }
        Column.push_back(E);
        //cout<<E<<", ";
    }
    cout<<endl;
     
    *ObjAux=model.get(GRB_DoubleAttr_ObjVal);
    //cout<<"Obj_Sol: "<<model.get(GRB_DoubleAttr_ObjVal)<<endl;
     
    for (int i=0; i<Nodes; i++) {
        delete[] Y[i];
    }
    delete[] Y;
     
    for (int i=0; i<Nodes+Targets; i++) {
        delete[] x[i];
    }
    delete[] x;
    delete env;
     
    return Column;
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
 
Vec Random_Sol(Mat &Adj,int Nodes, int T, int *Act_Sen, double alpha){
    //Declaration of solution vectos and required auxiliary variables
    Vec Solution;
    Vec Coverage;
    int RN, Cov_Lev=0;
     
    //initializaing vVectors
    for (int i=0; i<Nodes; i++) {
        Solution.push_back(0);
    }
     
    for (int i=0; i<T; i++) {
        Coverage.push_back(0);
    }
     
     
    //Generating random solution
     
    RND_DFS(Adj, Solution, Coverage, Nodes, T, 0 ,&Cov_Lev,  alpha);
    return Solution;
     
}
 
void RND_DFS(Mat &Adj, Vec &Sol,Vec &Cob, int Nodes, int T, int Cur_Node,int *CovTotal, double alpha ){
    int cont=1;
    int RDN_num;
     
    while (true) {
         
        Vec List_reach;
        for (int i=1; i<Nodes; i++) {
            if((Adj[Cur_Node][i]==1 or Adj[i][Cur_Node]==1) and (Sol[i]==0 or Sol[i]==1)) {
                List_reach.push_back(i);
            }
        }
        if (*CovTotal<T*alpha and List_reach.size()>0) {
            RDN_num=rand()%(List_reach.size());
            cont=List_reach[RDN_num];
            Sol[cont]=rand()%(1)+2;
            if (Sol[cont]==2) {
                for(int j=0;j<T;j++){
                    if(Adj[Nodes+j][cont]==1 and Cob[j]==0){
                        Cob[j]=1;
                        *CovTotal+=1;
                    }
                }
            }
            RND_DFS(Adj, Sol, Cob, Nodes, T, cont ,CovTotal,  alpha);
        }
        else{
            break;
        }
         
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
