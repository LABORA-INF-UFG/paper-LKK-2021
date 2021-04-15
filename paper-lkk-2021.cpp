#include <map>
#include <vector>
#include <ilcplex/ilocplex.h>
#include <fstream>
#include <iostream>
#include <cstdlib>  
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <time.h>
#include <limits.h>
#include <ctime>  

#include <cstring>
   
using namespace std;
   

  typedef struct {
  int source;
  int target;
  int data[1000][1000];
  map<pair<int, int>, IloBoolVar> variables;
   } Flow;
  int main (int argc, char **argv) {
  ////////////////////////////////////////////////////////////////
  double inicio = clock(); 
  ///////////////////////////////////////////////////////////////  
  // Usaremos FILE para a impressao 
  ofstream outputFile;
  outputFile.open("lpex1.txt");

  ////////////////////////////////////////////////////////////////////
  //Inicialização
  ///////////////////////////////////////////////////////////////////
  int instancia = 2;
  int numVertex = 100; //numero de nós 
  int numFlow = 80;  //quantidade de fluxos  
   vector<vector<int> > graph(numVertex + 1); 
   Flow *flows = new Flow[numFlow+1];
  int mo = 3; // quantidade inicial de nós e mo é igual também a quatidade de aresta acrescida em cada nó novo
  int t = numVertex - mo; //quantidade de nós que quero acrescentar
  int numEdge1= 2*mo*t; //quantidade de arestas do grafo
  int maxnumvar= 2*numEdge1*numFlow;//quantidade de variáveis
  int alpha0 = INT_MAX;
  int ksi = 3;//quantidade distinda de valores da qualidade   
  float *q = new float[ksi];//vetor que armazena os valores das qualidades distintas  em ordem crescente
   q[0]= 20;
   q[1]= 30;
   q[2]= 60;
   float theta1 = 1;
   float theta2 = 1;
   int itMax = numFlow*30;
  /////////////////////////////////////////////////////////
  //matriz de adjacência do grafo
  /////////////////////////////////////////////////////////
	int** inc=new int*[numVertex];
    for(int i = 0; i < numVertex; i++) {
      inc[i]=new int[numVertex];
     }
	
	
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Criar uma matriz X com (numFLow) linhas e (maxnumvar) colunas onde cada linha é o vetor solução (variables) de cada iteração
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   int** X = new int*[itMax];
   for (int i = 0;i < itMax;i++ ) {
       X[i] = new int[maxnumvar];
    }
  
    //////////////////////////////////////////////////////////////////////////////////////////////
  //Criar uma matriz  com (numVertex) linhas e (numVertex) colunas e (NumFlow) profundidade vetor solução (variables) de cada iteração
   /////////////////////////////////////////////////////////////////////////////////////////////
  float*** variableT = new float**[numVertex];//novo
  for (int i = 0; i < numVertex; ++i) {
    variableT[i] = new float*[numVertex];
	for (int j = 0; j < numVertex; ++j){
		variableT[i][j] = new float[numFlow];}
  }

   //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Criar uma matriz VO de ordem (numFLow x 2)  que em cada linha representa o vetor objetivo(z1, z2) candidatos a solução eficiente de P.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
   int** VO = new int*[itMax];
   for (int i = 0;i < itMax;i++ ) {
       VO[i] = new int[2];
   }
  
   //////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////Criar uma matriz PO de ordem (numFLow x 2)  que em cada linha representa o vetor objetivo(z1, z2)  solução eficiente de P.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
   int** PO = new int*[itMax];
   for (int i = 0;i < itMax;i++ ) {
       PO[i] = new int[2];
   }

   
   //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Matriz salto numFlow x numFlow que armazena o total de saltos de cada fluxo de cada iteração.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
   int** salto = new int*[itMax];
   for (int i = 0;i < itMax;i++ ) {
       salto[i] = new int[numFlow];
   }


   int** sumFlow=new int*[numVertex];// matriz dos pesos total de cada aresta
    for(int i = 0; i < numVertex; i++) {
      sumFlow[i]=new int[numVertex];
     }


   ///////////////////////////////////////////////////////////////////////////////////////////////
   int** loc = new int*[numVertex];
   for (int i = 0;i < numVertex;i++ ) {
       loc[i] = new int[numVertex];
   }
   //////////////////////////////////////////////////////////////////////////////////////////////
   float** prob = new float*[numVertex];
   for (int i = 0;i < numVertex;i++ ) {
       prob[i] = new float[numVertex];
   }

   //////////////////////////////////////////////////////////////////////////////////////////////
  //Determinar a matriz numVertex x numVertex x nunmFlow [s][d][f] os valores de wfsd
   /////////////////////////////////////////////////////////////////////////////////////////////
  float*** pesoTotal = new float**[numVertex];
  for (int i = 0; i < numVertex; ++i) {
    pesoTotal[i] = new float*[numVertex];
	for (int j = 0; j < numVertex; ++j){
		pesoTotal[i][j] = new float[numFlow];}
  }

  
   //////////////////////////////////////////////////////////////////////////////////////////////
  //Matriz que armazena a soma dos pesos dos fluxos que passaram na(s) aresta(s) gargalo 
    //////////////////////////////////////////////////////////////////////////////////////////////
   int** sumPesoBottleneck = new int*[itMax];
   for (int i = 0;i < itMax;i++ ) {
       sumPesoBottleneck[i] = new int[numEdge1];
    }
  //////////////////////////////////////////////////////////////////////////////////////////////
  //Matriz que armazena a soma dos pesos dos fluxos que passaram na(s) aresta(s) gargalo 
    //////////////////////////////////////////////////////////////////////////////////////////////
   int** sumQualiBottleneck = new int*[itMax];
   for (int i = 0;i < itMax;i++ ) {
       sumQualiBottleneck[i] = new int[numEdge1];
    }
   

   //////////////////////////////////////////////////////////////////////////////////////////////
  //Matriz que armazena a soma dos pesos dos fluxos que passaram na(s) em todas as arestas do grafo 
    //////////////////////////////////////////////////////////////////////////////////////////////
   int** sumPesoAll = new int*[itMax];
   for (int i = 0;i < itMax;i++ ) {
       sumPesoAll[i] = new int[numEdge1];
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
   //Matriz que armazena a soma dos pesos dos fluxos que passaram na(s) em todas as arestas do grafo 
    //////////////////////////////////////////////////////////////////////////////////////////////
   int** sumQualiAll = new int*[itMax];
   for (int i = 0;i < itMax;i++ ) {
       sumQualiAll[i] = new int[numEdge1];
    }


   /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  	int *variable0 = new int[maxnumvar]; //vetor solução do caminho mínimo com numvar elementos
    int *variable = new int[maxnumvar]; //vetor solução com numvar elementos
	int *bottleneck = new int[itMax]; // vetor com o gargalo em cada iteração após resolver o P-alpha em cada iteração********
	int *ori = new int[numFlow];//vetor que *contem a origem de cada fluxo
    int *dest = new int[numFlow]; //vetor que *contem o destino de cada fluxo
    float *Gamma = new float[numFlow]; //vetor que *contem o caminho mínimo de cada fluxo
    int  *peso = new int[numFlow]; //vetor que *contem os pesos dos fluxos
    int *sumPesoMax = new int[itMax];//vetor que armazena a soma maxima dos pesos dentre as aresta com gagalo naquela iteração
	int *sumPesoMaxAll = new int[itMax];//vetor que armazena a soma maxima dos pesos dos fluxos de todas arestas
	int *sumQualiMax = new int[itMax];//vetor que armazena a soma maxima dos valores da qualidade dentre as aresta com gagalo naquela iteração
	int *sumQualiMaxAll = new int[itMax];//vetor que armazena a soma maxima dos valores da qualidade dos fluxos de todas arestas

	///////////////////////////////////////////////////////////////////////////////////////
    //GRAFO
	/////////////////////////////////////////////////////////////////////////////////////
    int *conect = new int[numVertex]; //quantidade de arestas inseridas em cada nó
	int *nodeAdj = new int[numVertex];//vetor que *conta quantos nós adjacentes têm cada nó
	int sumConectTotal= 0;//total de conecção
	int *position = new int[numVertex+1]; //total de conecção até o nó i
	
	for(int l = 0; l < numVertex; ++l){
        nodeAdj[l] = 0;
     }

   
    for(int l = 0; l < numVertex; ++l){
        conect[l] = 0;
     }
	////////////////////////////////////////////////////////////////////////
  
    for (int i=0; i< numVertex; i++){
        for (int j = 0; j < numVertex;j++){
                  inc[i][j] = 0;
	    }
    }
	  
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int numEdge2 = 0;
    //srand(time(0)); //gera numeros aleatorios distintos
    srand(30);//mantém os mesmos numeros gerados aqui /////
	//ligar o nó mo+1 com (m) nós dos (mo) nós já existentes (aqui m = m0)
	
	for(int i = 1; i <= mo; i++){
		   graph[mo+1].push_back(i);
		   graph[i].push_back(mo+1);
		   inc[mo][i-1] = 1;
		   inc[i-1][mo] = 1;
	       nodeAdj[i-1] = nodeAdj[i-1] + 1;
		   nodeAdj[mo] = nodeAdj[mo];
	       conect[i-1] = conect[i-1] + 1;
		   conect[mo] = conect[mo] + 1;
		   numEdge2 = numEdge2 + 2;
	}
	int tempo=0;
	position[0] = 0;
	position[1]= conect[0];
	for(int i = 2; i <= mo+1; i++){
		position[i] = position[i-1]+conect[i-1];
	}
	sumConectTotal= position[mo+1];

	// calculando a probabilidade e a distancia minima entre as conectividades
    for(int i = 1; i <= mo+1; i++){
		prob[tempo][i-1]=  float(conect[i-1])/ (sumConectTotal);
	}
	
	//inserir novos nós

	int vmo; 

	/////////////////////////////////////////////////////////////////////////////////////////
	//inserindo o nó
	/////////////////////////////////////////////////////////////////////////////////////////
	int vpi = 0;
	 int novo = 0;
	for(int i = mo+2; i <= numVertex; i++){
		vmo = 0;
		while (vmo < mo){
			novo = (rand()%(sumConectTotal)) + 1;
			for(int j = 1; j <= i-1; ++j){
				if(inc[j-1][i-1] == 0 && inc[i-1][j-1] ==0 & novo <= position[j] && novo > position[j-1]){//j != vp[vmo-2]){
					graph[i].push_back(j);
					graph[j].push_back(i);
					vmo= vmo+1;
					conect[i-1]= conect[i-1]+1;
					conect[j-1]= conect[j-1]+1;
					inc[j-1][i-1] = 1;
		            inc[i-1][j-1] = 1; 
		            nodeAdj[i-1] = nodeAdj[i-1] + 1;
					nodeAdj[j-1] = nodeAdj[j-1] + 1;
					numEdge2=numEdge2 + 2;
				}
			}
		}
		//atualizando as posições 
		position[0] = 0;
	    position[1]= conect[0];
	    for(int j = 2; j <= i; j++){
		   position[j] = position[j-1]+conect[j-1];
	    }
		sumConectTotal = position[i];
		// calculando a probabilidade 
		tempo = tempo+1;
		for(int j = 1; j <= i; j++){
		    prob[tempo][j-1]= float (conect[j-1])/(sumConectTotal);
	    }
		
	}
   //Calcular o número máximo E MÍNIMO de nós adjacentes
  ///////////////////////////////////////////////////////////////////////
  int maxnodeAdj = 0;
  int posmaxEdge = 0;
  for (int l=0; l< numVertex; ++l){
            if (nodeAdj[l] > maxnodeAdj){ 
			   maxnodeAdj = nodeAdj[l];
			   posmaxEdge = l+1;}
   } 
  int posminEdge =0;
  int minnodeAdj = INT_MAX;
  for (int l=0; l< numVertex; ++l){
            if (nodeAdj[l] < minnodeAdj){ 
			   minnodeAdj = nodeAdj[l];
			   posminEdge = l+1;}
   }
     
   //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Matriz qualidade numVertex x numVertex que armazena a qualidade da aresta.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int  posq; 
  int** qualitfy = new int*[numVertex];
   for (int i = 0;i < numVertex;i++ ) {
       qualitfy[i] = new int[numVertex];
   }
   srand(30);
   for (int i = 0;i < numVertex;i++ ) {
	   for (int j = 0;j < numVertex;j++ ) {
		   if(inc[i][j] == 1){
		       posq = (rand()%(ksi)) + 1;
		       qualitfy[i][j] = q[posq-1];
			   qualitfy[j][i] = q[posq-1];
		   }
	       if(inc[i][j] == 0){
	            qualitfy[i][j] = INT_MAX;
				qualitfy[j][i] = INT_MAX;
		   }

		}
   }

	   
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   int contador;               // Declara a quantidade de instancias.
   contador=0;                 
      
   while (contador <= instancia)      
  {
      double inicio = clock();
	  contador++;  
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   int alpha = INT_MAX;
   int numEdge=0;
   int cont=0; //numero de vezes que P_alpha é viável
   int numdelete = 0; //quantidade de soluções dominadas
   int numPareto = 0; //quantidade de soluções  Pareto-ótimas
   // int prod = 0;
   //int pos = 0;
   //int cost = 0;
     
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   /*
    for (int i = 0; i < itMax; ++i){
        bottleneck[i] = 0;
    } 

    for (int i = 0; i < numEdge; ++i){
        sumPesoMax[i] = 0;
    }
   
	for(int i = 0; i < itMax; ++i){
	  for(int j = 0; j < numEdge; ++j){
	   sumPesoBottleneck[i][j] = 0; 
     }
  }
   for(int i = 0; i < itMax; ++i){
	  for(int j = 0; j < numEdge; ++j){
	   sumPesoAll[i][j] = 0; 
     }
   }
   
   for (int i = 0; i < numEdge; ++i){
        sumPesoMaxAll[i] = 0;
    }

	for (int i = 0; i < numEdge; ++i){
        sumQualiMax[i] = 0;
    }
   
   for(int i = 0; i < itMax; ++i){
	  for(int j = 0; j < numEdge; ++j){
	   sumQualiBottleneck[i][j] = 0; 
     }
  }
      
  for(int i = 0; i < itMax; ++i){
	  for(int j = 0; j < numEdge; ++j){
	   sumQualiAll[i][j] = 0; 
     }
  }
  for (int i = 0; i < numEdge; ++i){
        sumQualiMaxAll[i] = 0;
    }
  */
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Construção dos Fluxos dinamico
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////vetor que *contem a origem de cada fluxo///////////////////////////////////////////////
  for(int i = 0; i < numFlow; ++i){
	   ori[i] = 0; 
     }
  /////////////////////////vetor que *contem o destino de cada fluxo///////////////////////////////////////////////
  for(int i = 0; i < numFlow; ++i){
	     dest[i] = 0; 
       }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////mesmo destino e origens quaisquer ///////////////////////////////////////////////////////////////////////////////
  
  for(int i = 1;i <= numFlow;i++){
	       memset(flows[i].data, 0, sizeof(flows[i].data));
  }
    srand(time(0)); //gera numeros aleatorios distintos
   //srand(30);//mantém os mesmos numeros  gerados aqui /////
      int d = (rand()%(numVertex)) + 1;
	  for(int i = 1;i <= numFlow;i++){
	      ori[i]=d;
		  dest[i]=d;
		  while(ori[i] == dest[i]){
			ori[i] = (rand()%(numVertex)) + 1; 
		  }
		 flows[i].source = ori[i];
	     flows[i].target = dest[i];
	  }
  
    /////////////////////////////////////////////////////////
   ////////origens e destinos quaisquer///////////////////////////////////////////////////////////////////////////////
 
   /*for(int i = 1;i <= numFlow;i++){
	       memset(flows[i].data, 0, sizeof(flows[i].data));
  }
  
  srand(time(0)); //gera numeros aleatorios distintos
  for(int i = 1;i <= numFlow;i++){
	  while( ori[i-1] == dest[i-1] ){    
	       ori[i-1] = (rand()%(numVertex)) + 1; //gera um inteiro de 1 a numVertex
		   dest[i-1] = (rand()%(numVertex)) + 1; }
	   flows[i].source = ori[i-1];
	   flows[i].target = dest[i-1];}
  
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int i = 1;i <= numFlow;i++){
              cout <<"Fluxo"<< i << "- origem ="<< flows[i].source << " - destino =" << flows[i].target << endl;
	          //outputFile <<"Fluxo"<< i << "- origem ="<< flows[i].source << " - destino =" << flows[i].target << endl;
   }
 

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Determinar o caminho mínimo de cada fluxo
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //A construção do modelo
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  IloEnv   env;
  IloModel model(env);
  IloBoolVarArray variables(env);
  IloRangeArray constraints(env);
  IloNumArray vals(env);
  CPX_ON;//conserva a memória quando possível
  CPX_NODESEL_DFS; //busca em profundidade
   
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Construção das variáveis
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int numvar = 0;

  for (int i = 1; i <= numFlow; i++) {
    map<pair<int, int>, bool> added;
    for (int j = 1; j <= numVertex; j++) {
      for (int k = 0; k < graph[j].size(); k++) {
        stringstream varName;
        varName << "x_" << i << "_" << j << "_" << graph[j][k];
        IloBoolVar newVar = IloBoolVar(env, 0, 1, varName.str().c_str());
        flows[i].variables[make_pair(j, graph[j][k])] = newVar;
        variables.add(newVar);
		
		numvar=numvar+1; //quantidade de variáveis
      }
    }
  }
  
  numEdge=numvar/numFlow; //quantidade de arestas
       
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Criando as restrições . 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Esse loop a seguir é utilizado nas duas primeiras restriçoes
  for (int i = 1; i <= numFlow; i++) {
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //criando restrição 1: garante que cada fluxo saia da origem e que cada fluxo alcança seu destino
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	stringstream out0ConstraintName;
    stringstream in0ConstraintName;
    
    out0ConstraintName << "Out_flow_" << i;
    in0ConstraintName << "In_flow_" << i;

    IloExpr out0Constraint(env);
    IloExpr in0Constraint(env);

     for (int j = 0; j < graph[flows[i].source].size(); j++) {
          out0Constraint += flows[i].variables[make_pair(flows[i].source, graph[flows[i].source][j])];
          out0Constraint -=   flows[i].variables[make_pair(graph[flows[i].source][j], flows[i].source)];

    }

    for (int j = 0; j < graph[flows[i].target].size(); j++) {
         in0Constraint += flows[i].variables[make_pair(flows[i].target, graph[flows[i].target][j])];       
         in0Constraint -= flows[i].variables[make_pair(graph[flows[i].target][j], flows[i].target)];

    }
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//criando a restrição 2: garante que o tráfego de cada fluxo tenha um único caminho
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// esse loop abaixo exclui os vértices de origem e destino
    for (int j = 1; j <= numVertex; j++) {
      if (j == flows[i].source || j == flows[i].target)
        continue;

	  
      IloExpr path0Constraint(env);
      stringstream path0ConstraintName;

      path0ConstraintName << "Path_" << i << "_" << j;

      for (int k = 0; k < graph[j].size(); k++) {
        path0Constraint += flows[i].variables[make_pair(graph[j][k], j)];
        path0Constraint -= flows[i].variables[make_pair(j, graph[j][k])];
      }

      constraints.add(IloRange(env, 0, path0Constraint, 0, path0ConstraintName.str().c_str()));
       }

      constraints.add(IloRange(env, 1, out0Constraint, 1, out0ConstraintName.str().c_str()));
      constraints.add(IloRange(env, -1, in0Constraint, -1, in0ConstraintName.str().c_str()));
       }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Criando a restrição 3: garante a quantidade máxima de fluxos em cada aresta 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  for (int i = 1; i <= numVertex; i++) {
    for (int j = 0; j < graph[i].size(); j++) {//aquizero <
      
	  IloExpr alpha0Constraint(env);
      stringstream alpha0ConstraintName;

      alpha0ConstraintName << "Alpha_" << i << "_" << graph[i][j];
	  

      for (int k = 1; k <= numFlow; k++) {
        alpha0Constraint += flows[k].variables[make_pair(i, graph[i][j])];
		
      }


      constraints.add(IloRange(env, 0, alpha0Constraint, alpha0, alpha0ConstraintName.str().c_str()));

    }
  }

  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Criando o modelo para achar o caminho mínimo
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  model.add(variables);
  model.add(IloMinimize(env, IloSum(variables)));
  model.add(constraints);
  IloCplex cplex(model);
  cplex.setOut(env.getNullStream()); // Removendo logging do cplex
  cplex.solve();
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Imprimindo os resultados do modelo//////////////////////////////////////////////////////////////////////////////////////
   env.out() << "Solution status = " << cplex.getStatus() << endl;
   env.out() << "Solution value  = " << cplex.getObjValue() << endl;
    
 /* for(int i = 1; i <= numvar; ++i){
	     cplex.getValues(vals, variables);    
	    variable0[i-1] = cplex.getValue(variables[i-1]);
		//outputFile << "variable_ "<<i<<" = " << variable0[i-1] << endl;
   }
   */

 // Armazenar o total de saltos de cada fluxo (caminho mínimo) 
/*  for(int k = 0; k < numFlow; ++k){
		 for(int i = 0; i < numEdge; ++i){
		 Gamma[k] = Gamma[k] + variable0[i + k*numEdge]; 
     }
  }
  */
 for(int k = 1; k <= numFlow; ++k){
	 Gamma[k-1]=0;
 }

  for(int k = 1; k <= numFlow; ++k){
    for(int i = 1; i <= numVertex; ++i){
	   for(int j = 0; j < graph[i].size(); ++j){
			 cplex.getValues(vals, variables);
			 Gamma[k-1] = Gamma[k-1] + cplex.getValue(flows[k].variables[make_pair(i, graph[i][j])]);
		}
    }
   }


  //Limpar 
	  cplex.clearModel();
	  cplex.clear();
	  variables.end();
      constraints.endElements();
      constraints.end();
      cplex.end();
      model.end();
	  env.end();
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Achar caminho máximo
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  float gammaMax = 0;
  for(int k = 0; k < numFlow; ++k){
	         if(gammaMax < Gamma[k]){
				 gammaMax = Gamma[k];
			 }
  }
  //outputFile << "Gamma Maximo = "<< gammaMax << endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Determinar os pesos dos fluxos
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  float peso0;
  float auxpeso;
  for(int k = 0; k < numFlow; ++k){
	         peso0 = gammaMax/Gamma[k];
			 auxpeso = ceil(peso0);
			 peso[k] = auxpeso;
			 
  }
  
  //Achar o peso máximo
  float pesoMax = 0;
  for(int k = 0; k < numFlow; ++k){
	         if(pesoMax < peso[k]){
				 pesoMax = peso[k];
			 }
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Determinar os valores de betas
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  float betaTwo;
  float betaOne;
  float auxbeta = theta1;
  float auxq = q[ksi-1]/pesoMax;
  auxbeta= auxbeta*auxq;
  betaOne = ceil(auxbeta);
  betaTwo = theta2;
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Determinar os valores de W_f_s_d (peso do fluxo f atravessar a aresta e_sd)
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (int i = 0; i < numVertex; ++i) {
    	for (int j = 0; j < numVertex; ++j){
			if(inc[i][j]==1){
			   for(int k = 0; k < numFlow; ++k){
					pesoTotal[i][j][k] = betaOne*peso[k] + betaTwo*qualitfy[i][j];
			   }
			}
			if(inc[i][j]==0){
			   for(int k = 0; k < numFlow; ++k){
					pesoTotal[i][j][k] = INT_MAX;
			   }
			}
		}
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //Construção do Loop para a resolução de P de Alpha
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  do{     
 
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //A construção do modelo
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  IloEnv   env;
  IloModel model(env);
  IloBoolVarArray variables(env);
  IloRangeArray constraints(env);
  IloNumArray vals(env);
  CPX_ON;//conserva a memória quando possível
  CPX_NODESEL_DFS; //busca em profundidade
   
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Construção das variáveis
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int numvar = 0;

  for (int i = 1; i <= numFlow; i++) {
    map<pair<int, int>, bool> added;
    for (int j = 1; j <= numVertex; j++) {
      for (int k = 0; k < graph[j].size(); k++) {
        stringstream varName;
        varName << "x_" << i << "_" << j << "_" << graph[j][k];
        IloBoolVar newVar = IloBoolVar(env, 0, 1, varName.str().c_str());
        flows[i].variables[make_pair(j, graph[j][k])] = newVar;
        variables.add(newVar);
		
		numvar=numvar+1; //quantidade de variáveis
      }
    }
  }
  
  numEdge=numvar/numFlow; //quantidade de arestas
       
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Criando as restrições . 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Esse loop a seguir é utilizado nas duas primeiras restriçoes
  for (int i = 1; i <= numFlow; i++) {
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //criando restrição 1: garante que cada fluxo saia da origem e que cada fluxo alcança seu destino
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	stringstream outConstraintName;
    stringstream inConstraintName;
    
    outConstraintName << "Out_flow_" << i;
    inConstraintName << "In_flow_" << i;

    IloExpr outConstraint(env);
    IloExpr inConstraint(env);

     for (int j = 0; j < graph[flows[i].source].size(); j++) {
          outConstraint += flows[i].variables[make_pair(flows[i].source, graph[flows[i].source][j])];
          outConstraint -=   flows[i].variables[make_pair(graph[flows[i].source][j], flows[i].source)];

    }

    for (int j = 0; j < graph[flows[i].target].size(); j++) {
         inConstraint += flows[i].variables[make_pair(flows[i].target, graph[flows[i].target][j])];       
         inConstraint -= flows[i].variables[make_pair(graph[flows[i].target][j], flows[i].target)];

    }
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//criando a restrição 2: garante que o tráfego de cada fluxo tenha um único caminho
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// esse loop abaixo exclui os vértices de origem e destino
    for (int j = 1; j <= numVertex; j++) {
      if (j == flows[i].source || j == flows[i].target)
        continue;

	  
      IloExpr pathConstraint(env);
      stringstream pathConstraintName;

      pathConstraintName << "Path_" << i << "_" << j;

      for (int k = 0; k < graph[j].size(); k++) {
        pathConstraint += flows[i].variables[make_pair(graph[j][k], j)];
        pathConstraint -= flows[i].variables[make_pair(j, graph[j][k])];
      }

      constraints.add(IloRange(env, 0, pathConstraint, 0, pathConstraintName.str().c_str()));
       }

      constraints.add(IloRange(env, 1, outConstraint, 1, outConstraintName.str().c_str()));
      constraints.add(IloRange(env, -1, inConstraint, -1, inConstraintName.str().c_str()));
       }
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Criando a restrição 3: garante a quantidade máxima de fluxos em cada aresta considerando separado os fluxos q passam 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  for (int i = 1; i <= numVertex; i++) {
    for (int j = 0; j < graph[i].size(); j++) {//aquizero <
      
	  IloExpr alphaConstraint(env);
      stringstream alphaConstraintName;

      alphaConstraintName << "Alpha_" << i << "_" << graph[i][j];
	  

      for (int k = 1; k <= numFlow; k++) {
        alphaConstraint += pesoTotal[i-1][graph[i][j]-1][k-1]*flows[k].variables[make_pair(i, graph[i][j])];
		
      }


      constraints.add(IloRange(env, 0, alphaConstraint, alpha, alphaConstraintName.str().c_str()));

    }
  }
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Criando o modelo P_alpha
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  model.add(variables);
  model.add(IloMinimize(env, IloSum(variables)));
  model.add(constraints);
      
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //Imprimindo os dados e o modelo P_alpha
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////imprimindo os resultados
  if(cont == 0){
   cout<<"Number of nodes = "<< numVertex << endl;
   cout<<"Number of edges = "<< numEdge << "=" << numEdge2 << endl;
   cout<<"Number of flows = "<< numFlow << endl;
   cout<<"Number of variables = "<< numvar << endl;
  
    }
      
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Otimiza o problema P_alpha e obtém a solução
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      IloCplex cplex(model);
      cplex.setOut(env.getNullStream()); // Removendo logging do cplex
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   if ( cplex.solve() ) {
	   cont = cont + 1;}
   else {
	  env.error() << "Failed to optimize LP" << endl;
	  cout << "The P-Alpha is infeasible to iterations >=  " << cont+1 << endl;  
	  //outputFile << "The P-Alpha is infeasible to iterations >=  " << cont+1 << endl;
	  cout <<"The quantify of times that P_alpha was solved = " << cont+1 << endl;
      cout << "The quantity of solutions dominated = " << numdelete << endl;
      cout << "The quantity of optimal-Pareto solutions, i.e., |X| = " << numPareto << endl;
   	  double fim = clock();
      double  tempo = fim - inicio;
	  //outputFile << "The time of execution until here is  " << tempo / (double)CLOCKS_PER_SEC << " seconds" << endl;	
	  cout << "The time of execution until here is  " << tempo / (double)CLOCKS_PER_SEC << " seconds" << endl;
	  int primeiro = 0;
	  for(int i=0;i<=cont;i++){
		  if( PO[i][0] != INT_MAX && primeiro == 0 ){
			  cout << numPareto << "  " << numdelete << "  " << numPareto+ numdelete << "  " << cont+1  << "  " << tempo / (double)CLOCKS_PER_SEC << "  " << VO[0][0] << "  " << VO[cont-1][0] << "  " << sumPesoMax[0] << "  " << sumPesoMax[i] << "  " << sumPesoMax[cont-1] << "  " << PO[i][0] << "  " << PO[i][1]  <<  "  "  << PO[cont-1][0] << "  " << PO[cont-1][1]  <<endl; 
		      outputFile  << numPareto << "  " << numdelete << "  " << numPareto+ numdelete << "  " << cont+1  << "  " << tempo / (double)CLOCKS_PER_SEC << "  " << VO[0][0] << "  " << VO[cont-1][0] << "  " << sumPesoMax[0] << "  " << sumPesoMax[i] << "  " << sumPesoMax[cont-1] << "  " << PO[i][0] << "  " << PO[i][1]  <<  "  "  << PO[cont-1][0] << "  " << PO[cont-1][1]  <<endl; 
	          for (int k = 0; k < numFlow; k++ ){
				  outputFile <<salto[0][k] << "   " << salto[cont-1][k] << endl;}
			  for(int j=i; j<cont;j++){
				  if(PO[j][0] != INT_MAX){
					  if(j == i){
						  outputFile << "sumPesoMax" << "   " <<"sumPesoMaxAll" << "   " << "sumQualiMax" << "   " << "sumQualiMaxAll" << endl;}
					      outputFile << sumPesoMax[j] << "   " << sumPesoMaxAll[j] << "   " << sumQualiMax[j] << "   " << sumQualiMaxAll[j] << endl;}}
			  outputFile<<"*************"<< endl;
			  primeiro = primeiro + 1;

		  }
	  }
	  
	  ///Limpar 
	  cplex.clearModel();
	  cplex.clear();
	  variables.end();
      constraints.endElements();
      constraints.end();
      cplex.end();
      model.end();
	  env.end();
	  break;
   }
    
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Imprimindo os resultados do modelo//////////////////////////////////////////////////////////////////////////////////////
   env.out() << "Solution status = " << cplex.getStatus() << endl;
   env.out() << "Solution value  = " << cplex.getObjValue() << endl;
   
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Atualizando o valor de alpha
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int i = 1; i <= numVertex; ++i){
	   for(int j = 0; j < graph[i].size(); ++j){
		 for(int k = 1; k <= numFlow; ++k){
			 cplex.getValues(vals, variables);
			 variableT[i-1][graph[i][j]-1][k-1] = cplex.getValue(flows[k].variables[make_pair(i, graph[i][j])]);
		}
    }
   } 


   for(int i = 0; i < numvar; ++i){
	    variable[i] = cplex.getValue(variables[i]);
   }
   for(int i = 0; i < numVertex; ++i){
	   for(int j = 0; j < numVertex; ++j){
		      sumFlow[i][j] = 0;
        }
   }
   
    for(int i = 0; i < numVertex; ++i){
	   for(int j = 0; j < numVertex; ++j){
		   if(inc[i][j]==1){
		       for(int k = 0; k < numFlow; ++k){
			  	       sumFlow[i][j] += pesoTotal[i][j][k]*variableT[i][j][k];
		       }
		   }
       }
   }
     
   bottleneck[cont-1] = 0;
for(int i = 0; i < numVertex; ++i){
	   for(int j = 0; j < numVertex; ++j){
		   if(inc[i][j]==1 && sumFlow[i][j] > bottleneck[cont-1]){ 
				 bottleneck[cont-1] = sumFlow[i][j];
		   }
    } 
 }

  alpha = bottleneck[cont-1] - 1;  //o novo valor de alpha
  
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Armazenar o total de saltos de cada fluxo da iteração atual 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 for(int k = 0; k < numFlow; ++k){
	salto[cont-1][k] = 0;}

  for(int k = 0; k < numFlow; ++k){
    for(int i = 0; i < numVertex; ++i){
	   for(int j = 0; j < numVertex; ++j){
			 if(inc[i][j]==1){
		      salto[cont-1][k] = salto[cont-1][k] + variableT[i][j][k];
			 }
		}
    }
   }
 
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Armazenar o total de pesos da aresta(s) gargalo da iteração atual 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int contEdge = 0;
  for(int i = 0; i < numVertex; ++i){
	 for(int j = 0; j < numVertex; ++j){
		 if(inc[i][j] == 1 && sumFlow[i][j] == bottleneck[cont-1]){
			 sumPesoBottleneck[cont-1][contEdge] = 0; 
			 for(int k = 0; k < numFlow; ++k){
					sumPesoBottleneck[cont-1][contEdge] = sumPesoBottleneck[cont-1][contEdge] + peso[k]*variableT[i][j][k];
			}
		    //outputFile<<"soma dos pesos_"<<cont-1 <<"_"<<contEdge << " = " << sumPesoBottleneck[cont-1][contEdge]<<endl;
		    contEdge= contEdge +1;

		}
     }
   }

  sumPesoMax[cont-1]=0;
  for(int k = 0; k < contEdge; ++k){
	if(sumPesoBottleneck[cont-1][k] > sumPesoMax[cont-1]){
		sumPesoMax[cont-1] = sumPesoBottleneck[cont-1][k];
	}
  }
  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int contEdge2 = 0;

  
  for(int i = 0; i < numVertex; ++i){
	 for(int j = 0; j < numVertex; ++j){
		 if(inc[i][j] == 1 ){
			 sumPesoAll[cont-1][contEdge2] = 0; 
			 for(int k = 0; k < numFlow; ++k){
					sumPesoAll[cont-1][contEdge2] = sumPesoAll[cont-1][contEdge2] + peso[k]*variableT[i][j][k];
		     }
		 
		 // outputFile<<"soma dos pesosT_"<<cont-1 <<"_"<<contEdge2 << " = " << sumPesoAll[cont-1][contEdge2]<<endl;
		  contEdge2= contEdge2 +1;
		 }
		}
   }

  sumPesoMaxAll[cont-1]=0;
  for(int k = 0; k < numEdge; ++k){
	if(sumPesoAll[cont-1][k] > sumPesoMaxAll[cont-1]){
		sumPesoMaxAll[cont-1] = sumPesoAll[cont-1][k];
	}
  }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Armazenar o total dos valores das qualdidades da aresta(s) gargalo da iteração atual 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int contEdgeQ = 0;
  
  for(int i = 0; i < numVertex; ++i){
	 for(int j = 0; j < numVertex; ++j){
		 if(inc[i][j] == 1 && sumFlow[i][j] == bottleneck[cont-1]){
			 sumQualiBottleneck[cont-1][contEdgeQ] = 0;
			 for(int k = 0; k < numFlow; ++k){
					sumQualiBottleneck[cont-1][contEdgeQ] = sumQualiBottleneck[cont-1][contEdgeQ] + qualitfy[i][j]*variableT[i][j][k];
			  }
		     //outputFile<<"Sum-Quali " << qualitfy[i-1][graph[i][j]-1] << "_" <<cont-1 <<"_"<<contEdgeQ << " = " << sumQualiBottleneck[cont-1][contEdgeQ]<<endl;
		     contEdgeQ= contEdgeQ +1;

		}
     }
   }

  sumQualiMax[cont-1]=0;
  for(int k = 0; k < contEdgeQ; ++k){
	if(sumQualiBottleneck[cont-1][k] > sumQualiMax[cont-1]){
		sumQualiMax[cont-1] = sumQualiBottleneck[cont-1][k];
	}
  }
  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int contEdgeQ2 = 0;
  for(int i = 0; i < numVertex; ++i){
	 for(int j = 0; j < numVertex; ++j){
		 if(inc[i][j] == 1){
		 	 sumQualiAll[cont-1][contEdgeQ2] = 0; 
			 for(int k = 0; k < numFlow; ++k){
				sumQualiAll[cont-1][contEdgeQ2] = sumQualiAll[cont-1][contEdgeQ2] + qualitfy[i][j]*variableT[i][j][k];

			}
		 //outputFile<<"soma dos valores das qualidadesT_"<< cont-1 <<"_"<<contEdgeQ2 << " = " << sumQualiAll[cont-1][contEdgeQ2]<<endl;
		  contEdgeQ2= contEdgeQ2 +1;
		 }
		}
   }

  sumQualiMaxAll[cont-1]=0;
  for(int k = 0; k < numEdge; ++k){
	if(sumQualiAll[cont-1][k] > sumQualiMaxAll[cont-1]){
		sumQualiMaxAll[cont-1] = sumQualiAll[cont-1][k];
	}
  }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Armazenar a solução (variables) da iteração atual  na linha cont-1 da matriz X com (numFLow)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 for(int j = 0; j < numvar; ++j){
        X[cont-1][j] = variable[j]; 
   }


  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Armazenar o vetor objetivo da iteração cont na matriz VO de ordem  numFlow x 2 na linha (cont-1) que representa em cada linha o vetor objetivo (f1,f2).
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    VO[cont-1][0] = bottleneck[cont-1];
    VO[cont-1][1] = cplex.getObjValue();
	VO[cont][0] = INT_MAX;
    VO[cont][1] = INT_MAX;

	PO[cont-1][0] = VO[cont-1][0];
    PO[cont-1][1] = VO[cont-1][1];
	PO[cont][0] = INT_MAX;
    PO[cont][1] = INT_MAX;
    
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Cleanup
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//freeGAMSModel();//preserves memory by dumping the GAMS model instance representation temporarily to disk
	cplex.clearModel();
	cplex.clear();
	cplex.clearModel();
	variables.end();
    constraints.endElements();
    constraints.end();
    cplex.end();
    model.end();
	env.end();
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "The objective vector [f1(x)  f2(x)] in the iteration " << cont << " =  [" << VO[cont-1][0] << "  " << VO[cont-1][1]  << "]" <<  endl; 
    //outputFile << "The objective vector [f1(x)  f2(x)] in the iteration " << cont << " =  [" << VO[cont-1][0] << "  " << VO[cont-1][1]  << "]" <<  endl; 
    cout << "The value of alpha in the next iteration , ie , iteration " << cont+1 << " = " << alpha << endl;
    //outputFile << "The value of alpha in the next iteration , ie , iteration  " << cont+1 << "= " << alpha << endl;
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Retirando as soluções dominada de X 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // quando  solução da iteração i domina a solução da iteração i-1 a linha i-1 de X é substituída por outra
  // com elementos nulos .  
   for ( int i = 0; i < numvar; i++ ){
         if( cont > 1 && VO[cont-1][1] == VO[cont-2][1])
             X[cont-2][i] = 0;
   }
 
  // quando  solução da iteração i domina a solução da iteração i-1 a linha i-1 de VO é substituída por outra
  // com elementos infinitos.

   for ( int i = 0; i < 2; i++ ){
	   if( cont > 1 && VO[cont-1][1] == VO[cont-2][1]){
             PO[cont-2][i] = INT_MAX;
			 
	   }
   }
   if( cont > 1 && VO[cont-1][1] == VO[cont-2][1]){
	   numdelete = numdelete+1; //calcula quantas soluções foram dominadas, ou seja, deletadas}
    }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   for (int i = 0; i <= cont; i++ ){
	 cout <<  "VO_" << i << " * " << VO[i][0] << "  " << VO[i][1] << endl;
     //outputFile <<  "VO_" << i << " * " << VO[i][0] << "  " << VO[i][1]  << endl;
    }
	   
	//imprime a matriz com os vetores Pareto ótimos
    cout << " The Matrix PO presents in each row the objective vector of Pareto-opitmal solutions, i.e, " << endl;
     //outputFile << " The Matrix PO presents in each row the objective vector of Pareto-opitmal solutions, i.e, " << endl;
   
   for (int i = 0; i <= cont; i++ ){
	 cout <<  "PO_" << i << " * " << PO[i][0] << "  " << PO[i][1] << endl;
     //outputFile <<  "PO_" << i << " * " << PO[i][0] << "  " << PO[i][1]  << endl;
   }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //Calcular e Imprimir no final o numero de soluções excluidas, de soluções não dominadas e o numero de iterações resolvidas.
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    numPareto = cont - numdelete; //quantidade de soluções Pareto ótimas
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   if (alpha == 0){// betaOne*pesoMax + betaTwo*q[ksi] ){
		cout << "The P-Alpha reached the last iteration = " << cont << endl;///atingiu o valor minimo de alpha em cont
		//outputFile << "The P-Alpha reached the last iteration = " << cont << endl;
        cout <<"The quantify of times that P_alpha was solved = " << cont << endl;
        cout << "The quantity of solutions dominated = " << numdelete << endl;
        cout << "The quantity of optimal-Pareto solutions, i.e., |X| = " << numPareto << endl;
        double fim = clock();
        double  tempo = fim - inicio;
	    int primeiro = 0;
	    for(int i=0;i<cont;i++){//<=cont
		  if( PO[i][0] != INT_MAX && primeiro == 0 ){
			  //cout << numPareto << "  " << numdelete << "  " << numPareto+ numdelete << "  " << cont+1  << "  " << tempo / (double)CLOCKS_PER_SEC << "  " << VO[0][0] << "  " << VO[cont-1][0] << "  " <<PO[i][0] << "  " << PO[i][1]  <<  "  "  << PO[cont-1][0] << "  " << PO[cont-1][1]  <<endl; 
		      cout << numPareto << "  " << numdelete << "  " << numPareto+ numdelete << "  " << cont+1  << "  " << tempo / (double)CLOCKS_PER_SEC << "  " << VO[0][0] << "  " << VO[cont-1][0] << "  " << sumPesoMax[0] << "  " << sumPesoMax[i] << "  " << sumPesoMax[cont-1] << "  " << PO[i][0] << "  " << PO[i][1]  <<  "  "  << PO[cont-1][0] << "  " << PO[cont-1][1]  <<endl; 
			  //outputFile<< numPareto << "  " << numdelete << "  " << numPareto+ numdelete << "  " << cont+1  << "  " << tempo / (double)CLOCKS_PER_SEC << "  " << VO[0][0] << "  " << VO[cont-1][0] << "  " <<PO[i][0] << "  " << PO[i][1]  <<  "  "  << PO[cont-1][0] << "  " << PO[cont-1][1]  <<endl; 
			  outputFile << numPareto << "  " << numdelete << "  " << numPareto+ numdelete << "  " << cont+1  << "  " << tempo / (double)CLOCKS_PER_SEC << "  " << VO[0][0] << "  " << VO[cont-1][0] << "  " << sumPesoMax[0] << "  " << sumPesoMax[i] << "  " << sumPesoMax[cont-1] << "  " << PO[i][0] << "  " << PO[i][1]  <<  "  "  << PO[cont-1][0] << "  " << PO[cont-1][1]  <<endl; 
			  
			  for (int k = 0; k < numFlow; k++ ){
				  outputFile <<salto[0][k] << "   " << salto[cont-1][k] << endl;}
			  for(int j=i; j<cont;j++){//novo
				  if(PO[j][0] != INT_MAX){
					  		  if(j == i){
						          outputFile << "sumPesoMax" << "   " <<"sumPesoMaxAll" << "   " << "sumQualiMax" << "   " << "sumQualiMaxAll" << endl;}
					          outputFile << sumPesoMax[j] << "   " << sumPesoMaxAll[j] << "   " << sumQualiMax[j] << "  " << sumQualiMaxAll[j] << endl;}}
			  outputFile<<"*************"<< endl;
			  primeiro = primeiro + 1;
		  }
	  }
   }

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
     double fim = clock();
     double  tempo = fim - inicio;
	 //outputFile << "The running time until here is  " << tempo / (double)CLOCKS_PER_SEC << " seconds" << endl;	
	 cout << "The running time  until here is  " << tempo / (double)CLOCKS_PER_SEC << " seconds" << endl;
	 
	 
	
	} while (alpha != 0);
  

  } 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "The end." << endl;   
   system ("pause");
   outputFile.close();
      
     delete [] nodeAdj;
     delete [] ori;
     delete [] dest;
     delete [] flows;
	 delete [] variable;
     delete [] bottleneck;
	 delete [] conect;
	 delete [] position;
	 delete [] Gamma;
	 delete [] variable0;
	 delete [] q;
	 delete [] peso;
	 delete [] sumPesoMax;
     delete [] sumPesoMaxAll;
     delete [] sumQualiMax;
     delete [] sumQualiMaxAll;
	 
	for(int i=0 ; i < itMax;i++) {
      delete [] X[i];}
      delete [] X;
 
    for(int i=0 ; i < itMax;i++) {
       delete [] PO[i];}
       delete [] PO;
   
	 
	for(int i=0 ; i < itMax;i++) {
     delete [] VO[i];}
     delete [] VO;
	 
	 
	 for(int i=0 ; i < numVertex;i++) {
      delete [] inc[i];}
      delete [] inc;

     for(int i=0 ; i < numVertex;i++) {
      delete [] loc[i];}
      delete [] loc;
	 
	 for(int i=0 ; i < numVertex;i++) {
      delete [] prob[i];}
      delete [] prob;

    for(int i=0 ; i < itMax;i++) {
     delete [] salto[i];}
     delete [] salto;
  
	for(int i=0 ; i < numVertex;i++) {
      delete [] qualitfy[i];}
      delete [] qualitfy;

	for(int i=0 ; i < numVertex;i++) {
      delete [] sumFlow[i];}
      delete [] sumFlow;
	
    for (int i = 0; i < numVertex; ++i) {
      for (int j = 0; j < numVertex; ++j){
		  delete [] pesoTotal[i][j];}
	  delete [] pesoTotal[i];}
    delete [] pesoTotal;

  for (int i = 0; i < numVertex; ++i) {
      for (int j = 0; j < numVertex; ++j){
		  delete [] variableT[i][j];}
	  delete [] variableT[i];}
    delete [] variableT;


	for (int i = 0;i < itMax;i++ ) {
		delete [] sumPesoBottleneck[i];}
	delete [] sumPesoBottleneck;
    
	for (int i = 0;i < itMax;i++ ) {
		delete [] sumPesoAll[i];}
	delete [] sumPesoAll;
    
	for (int i = 0;i < itMax;i++ ) {
		delete [] sumQualiBottleneck[i];}
	delete [] sumQualiBottleneck;
    
	for (int i = 0;i < itMax;i++ ) {
		delete [] sumQualiAll[i];}
	delete [] sumQualiAll;

   return 0;
} 