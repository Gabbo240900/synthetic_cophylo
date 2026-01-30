#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>

/*_________________________________________________________*/
/**
 * @author Nicolas Alcala;  
 *
 */

typedef int bool;
#define true 1
#define false 0

void shuffle(int *array, int n)
{
    if (n > 1) 
    {
        int i;
        for (i = 0; i < n - 1; i++) 
        {
          int j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

float rexp(float lambda){
  if(lambda==0) return(100000000); // a very large number (infinity)
  float R;
  R = (float)rand()/(double)(RAND_MAX); //random number between 0 and 1.
  return  -log(R)/lambda;               //transform into exponential variable
}

int rdist(float * dist, int nd){ //generate a random int from a specified distribution
  float p = rand()/(float)RAND_MAX; //determines how many H infected by 1st P
  float cumuld=0;
  int j;
  for(j=0;j<nd;j++){
    cumuld += dist[j];
    if(cumuld>p){
      return(j+1);
    }
  }
  return(nd);
}

int bd_P_cospe(int* H_edge, float* H_edge_length, int** pedges, int* nedges, bool* edge_alive, int NH, float l, float m, float p, int nd, float* d, float ls, float pc, float Tmrca,int* P_edge, float* P_edge_length,int N){
  /*#################### parameters #############################
  #H_edge: host edge list
  #H_edge_length: host edge length list
  #pedges: list of hosts associated with each parasite
  #nedges: number of hosts associated with each parasite
  #edge_alive: is each parasite species alive or extinct
  #NH: number of host species
  #l: parasite speciation rate
  #m: parasite extinction rate
  #p: probability of vertical transmission of parasite to both children hosts
         #- 0: parasite transmitted to single child host (specialist) 
         #- 1: parasite transmitted to both children hosts (generalist)
  #nd: maximum number of hosts per parasite
  #d: distribution of the number of hosts per "new" parasite lineage
  #ls: rate of Host switch
  #pc: cospeciation probability, given host speciation
  #TMRCA: parasite TMRCA
  #P_edge: parasite edges
  #P_edge_length: parasite edge length
  #N: number of tip parasites
  #############################################################*/
  int i,j,k,ll,o,il,jl,ih;
  for(j=0;j<2*N;j++){ P_edge_length[j] = 0;}
  for(j=0;j<4*N;j++){ P_edge[j] = 0;}
	
  //Build a list of all speciation events occurring on the Host phylogeny with their timing
  float Ths[NH-1]; //timing of all H spe events
  int** nlistt = malloc( (NH-1)*sizeof(int*)); // new int*[NH-1]; //list of H edges after all events
  for(j=0;j<NH-1;j++) nlistt[j] = malloc( (j+2)*sizeof(int));//new int[j+2]; //list of H edges after event j
  float** tcur = malloc( (NH-1)*sizeof(float*));//new float*[NH-1]; //list of H edge lengths after all events
  for(j=0;j<NH-1;j++) tcur[j] = malloc( (j+2)*sizeof(float)); //new float[j+2]; //list of H edge lengths after event j

  int children=0;
  
  for( j=0;j<NH*2;j++){
    if( H_edge[j*2] == H_edge[0]){
      if(children==0){
	nlistt[0][0]  = j; //add new lineages
	tcur[0][0]    = H_edge_length[j]; //add length of new lineages
      }else{
	nlistt[0][1] = j; //add new lineages
	tcur[0][1]   = H_edge_length[j]; //add length of new lineages
      }
      children++;
      if(children==2) break;
    }
  }  
  int edgeevlist[NH-2]; //list of edges where event occured
  for( i=0; i<NH-2; i++) edgeevlist[i] = -1;
  
  for( i=0; i<NH-2; i++){
    int wh=0;
    for( j=1;j<i+2;j++){
      if(tcur[i][wh]>tcur[i][j] ) wh = j; //find lineage where next event occurs (shortest branch)
    }
    Ths[i]        = tcur[i][wh]; //add timing of event
    edgeevlist[i] = nlistt[i][wh]; //edge with event
    for( j=0;j<i+2;j++){
      tcur[i+1][j]   = tcur[i][j] - tcur[i][wh]; //reduce length of edges where no event occurred
      nlistt[i+1][j] = nlistt[i][j]; //duplicate edges id
    }
    children=0;
    for( j=0;j<NH*2;j++){
      if( H_edge[j*2] == H_edge[2*nlistt[i][wh]+1]){
	if(children==0){
	  nlistt[i+1][wh]  = j; //add new lineages
	  tcur[i+1][wh]    = H_edge_length[j]; //add length of new lineages
	}else{
	  nlistt[i+1][i+2] = j; //add new lineages
	  tcur[i+1][i+2]   = H_edge_length[j]; //add length of new lineages
	}
	children++;
	if(children==2) break;
      }
    }  
  }
  Ths[NH-2] = tcur[NH-2][0]; //add length of last edge
  
  int pstart = 0; //which event occurred right after the apparition of 1st Parasite
  double Thtmp = 0; //current time
  for( i=NH-2; i>=0; i--){
    Thtmp += Ths[i];
    if(Thtmp>Tmrca){
      pstart = i;
      break;
    }
  }
    
  int P_tip_label[N]; //id of P tips
  int Nnode; //number of nodes

  for( j=0;j<2*N;j++){
    for( k=0;k<NH;k++){
      pedges[j][k] = 0;
    }
  }
  for( j=0;j<2*N;j++) nedges[j] = 0;
  
  nedges[0] = rdist(d,nd);  //number of H infected by 1st P
  if(nedges[0]>pstart+2) nedges[0] = pstart+2;
  printf("\tnumber of H infected by 1st P: %d\n",nedges[0]);

  //draw without replacement nh0 H
  int drawtmp[pstart+2]; //possible edges
  for (i=0; i<pstart+2; ++i){
    drawtmp[i] = nlistt[pstart][i]; //set values
  }
  shuffle(drawtmp ,pstart+2);//&drawtmp[0],&drawtmp[pstart+2]); // using built-in random generator:

  printf("\tfirst infected H: ");
  for( j=0; j<nedges[0];j++){
    pedges[0][j] = drawtmp[j]; //change
    printf("H%d, ",H_edge[2*pedges[0][j]+1]);
  }						
  printf("\n");

  int ed  = 0;     //current node ID
  int pt[N];    //current parasites
  for( j=1;j<N;j++) pt[j] = -1;
  pt[0]   = ed;
  int is   = 0;
  int npt = 1;
  int extP;        //extinct parasites
  
  for( j=0;j<N*2;j++) edge_alive[j]=false;
  float tp,tpp=0;
 
  bool is_tip[N*2];
  for( j=0;j<2*N;j++) is_tip[j]=false;
  int ntip = 0;
  
  /****** P phylogeny simulation loop  *****/
  float th; //time until next H speciation
  float thh=0,tppp=0;
  int npopS = nedges[0];
  for( ih=pstart; ih<NH-1; ih++){
    if((npt==0)||(2*is>=N)){
      return 0;
    }
    if(ih==pstart){ 
      th = Ths[ih] - (Thtmp- Tmrca); //time to next H speciation event
    }else{ 
      th = Ths[ih];
    }
    thh += th;
    
    if(ih==pstart){
      tp  = 0;
      tpp = 0;
    }else{
      tp   = rexp(npt*(l+m)+npopS*ls); 
      tpp  = tp;
    }

    while(tpp<th){
      if((npt==0)||(2*is>=N)) return 0;
      
      float event;
      if( (is==0) ){ event = 0;}
      else{ event = (float)(rand()*(npt*(l+m)+npopS*ls)/(double)(RAND_MAX)); }//which type of event: 1. spe 2. ext 3. H switch
      if(event<=npt*l){  //speciation
	if(2*(is+1)>=N) return(0);
	int pts = rand()%npt; //choose P which speciates
	//add children P to edge list
	P_edge[is*4]   = pt[pts]; //change size
	P_edge[is*4+1] = ed+1; //change size
	P_edge[is*4+2] = pt[pts]; //change size
	P_edge[is*4+3] = ed+2; //change size
	//update P edge length
	for( j=0;j<2*N;j++){ if(edge_alive[j]){P_edge_length[j] += tp;} }
	tppp   += tp;
	//update extinction status
	if(is>0) edge_alive[pt[pts]-1] = false; //speciating P is considered extinct
	edge_alive[ed] = true; //new P become alive
	edge_alive[ed+1] = true; //new P become alive
	for( j=0;j<nedges[pt[pts]];j++){
	  pedges[ed+1][j] = pedges[pt[pts]][j]; //"old" P keeps same H as parent
	}
	nedges[ed+1] = nedges[pt[pts]];
	int nhednew = rdist(d,nd); //random number of H for the new P
	if(nhednew >= nedges[pt[pts]]){//child P is transmitted to all H of parent
	  nedges[ed+2] = nedges[pt[pts]];
	  npopS += nedges[ed+2];
	  for( j=0;j<nedges[pt[pts]];j++){
	    pedges[ed+2][j] = pedges[pt[pts]][j];
	  }
	}else{//child P is transmitted to a subset of the H of parent
	  nedges[ed+2] = nhednew;
	  npopS += nedges[ed+2];
	  //draw without replacement nedges[ed+2] H
	  int drawtmp[nedges[pt[pts]]]; //possible edges
	  for (i=0; i<nedges[pt[pts]]; ++i) drawtmp[i]=pedges[pt[pts]][i]; //set values
	  shuffle(drawtmp,nedges[pt[pts]]);//random_shuffle( &drawtmp[0],&drawtmp[nedges[pt[pts]]]); //using built-in array shuffling
	  for( j=0; j<nedges[ed+2];j++){
	    pedges[ed+2][j] = drawtmp[j]; 
	  }	
	}//end else

	//update list of P
	pt[pts] = ed+1;//replace id
	pt[npt] = ed+2;//add id
	ed      = ed+2;//new last P id
	npt++;//increase number of parasites
	is++;//increase number of speciation events
      }else{
	if(event<npt*(l+m)){//extinction
	  for( j=0;j<2*N;j++){ if(edge_alive[j]){P_edge_length[j] += tp;} }
	  tppp   += tp;
	  int pte = rand()%npt;
	  edge_alive[pt[pte]-1] = false;	  
	  is_tip[pt[pte]] = true;
	  pt[pte] = pt[npt-1]; 	  
	  npt--;
	  if((npt==0)||(2*is>=N)) return 0;
	  ntip++;
	  npopS -= nedges[pt[pte]];//the extinct P cannot switch H
	}else{
	  if(event<=(npt*(l+m)+npopS*ls)){//H switch
	    for( j=0;j<2*N;j++){ if(edge_alive[j]){P_edge_length[j] += tp;} }
	    tppp   += tp;
	    float pP = npopS*(float)(rand()/(double)(RAND_MAX)); //to determine which P switches
	    int fromP = -1, toH = -1; //P which switches
	    int drawH[ih+2]; //possible edges
	    int ncum=0;
	    for( k=0;k<npt;k++){
	      ncum += nedges[pt[k]];
	      if( pP<=ncum ){//chosen P
		fromP = k;
		if(nedges[pt[k]]<ih+2){//HS only if P does not already infect all H
                  for (ll = 0; ll < ih + 2; ++ll) {
                    drawH[ll] = nlistt[ih][ll]; // set values
                  }
                  shuffle(drawH, ih + 2); // randomize candidate hosts

                  // pick the first host not already infected by this parasite
                  toH = -1;
                  for (int cand = 0; cand < ih + 2; cand++) {
                    int candidate = drawH[cand];
                    bool already = false;
                    for (o = 0; o < nedges[pt[fromP]]; o++) {
                      if (candidate == pedges[pt[fromP]][o]) {
                        already = true;
                        break;
                      }
                    }
                    if (!already) {
                      toH = candidate;
                      break;
                    }
                  }

                  // If every host is already infected (should be rare given the earlier check), skip this switch.
                  if (toH >= 0) {
                    pedges[pt[fromP]][nedges[pt[fromP]]] = toH;
                    nedges[pt[fromP]]++;
                    npopS++;
                  }
                  break;
		}// if P does not infect every H
	      }//chosen P
	    }
	  }//HS
	}
      }
      tp     = rexp(npt*(l+m)+npopS*ls); 
      tpp    += tp;             //total time
    }
    
    // H speciation event
    if((npt>0)&(ih<NH-2) ){ 
      int newH1=-1,newH2=-1; //new H
      for( k=0; k<NH*2; k++){
	if( (H_edge[2*k]==H_edge[2*edgeevlist[ih]+1])&(newH1==-1) ) newH1 = k;
	else{if( (H_edge[2*k]==H_edge[2*edgeevlist[ih]+1])&(newH1>-1) ) newH2 = k;}
      }
      for( il=0; il<npt; ++il){ //for each surviving P
	for( jl=0; jl<nedges[pt[il]]; ++jl){//for each infected H
	  if(edgeevlist[ih]==pedges[pt[il]][jl]){//if P infected speciating H
	    float rco = (rand()/(float)RAND_MAX); //cospeciation or not
	    if(rco>pc){// no cospeciation
	      float type = rand()/(float)RAND_MAX; //random type of P transmission after H speciation
	      if(type<d[0]){ //randomly attribute P to a single child H => specialist P
		if( (rand()/(float)RAND_MAX)<0.5 ) pedges[pt[il]][jl] = newH1;
		else pedges[pt[il]][jl] = newH2;
	      }else{
		pedges[pt[il]][jl]             = newH1;
		pedges[pt[il]][nedges[pt[il]]] = newH2;
		nedges[pt[il]]++;
		npopS++;
	      }//type 1
	      break; //each P can only infect a H once 
	    }else{//cospeciation
	      if(2*(is+1)>=N) return(0);
	      //add children P to edge list
	      P_edge[is*4]   = pt[il]; //change size
	      P_edge[is*4+1] = ed+1; //change size
	      P_edge[is*4+2] = pt[il]; //change size
	      P_edge[is*4+3] = ed+2; //change size

	      P_edge_length[pt[il]-1] += th - tpp+tp; //update branch lengths
	      P_edge_length[ed] -= th - tpp+tp; //update branch lengths
	      P_edge_length[ed+1] -= th - tpp+tp; //update branch lengths
	      	      
	      //update extinction status  // <- maybe recheck here if not update too soon
	      edge_alive[pt[il]-1] = false; //speciating P is considered extinct
	      edge_alive[ed] = true; //new P become alive
	      edge_alive[ed+1] = true; //new P become alive
	      
	      for( j=0;j<nedges[pt[il]];j++){
		if(j!=jl){
		  pedges[ed+1][j] = pedges[pt[il]][j]; //"old" P keeps same H as parent
		}else{
		  pedges[ed+1][j] = newH1; //add new H
		}
	      }
	      pedges[ed+2][0] = newH2;
	      nedges[ed+1] = nedges[pt[il]];
	      nedges[ed+2] = 1;
	      npopS++;
	      
	      pt[il] = ed+1;//replace id
	      pt[npt] = ed+2;//add id
	      
	      ed = ed+2; //increment assigned node ID
	      npt++;//increase number of parasites
	      is++;//increase number of speciation events
	    }//cospe  
	  }//type selection
	}// if P on speciating H
      }//for each surviving P
    }//if P
    
    for( j=0;j<2*N;j++){ if(edge_alive[j]){P_edge_length[j] += th - tpp+tp;}} //update branch lengths
    tppp   += th - tpp+tp; //update total time
  }
  
  for( j=0;j<N*2;j++){ if(edge_alive[j]){is_tip[j+1] = true;ntip++;}}
 
  //reordering the nodes for R phylo class
  int newnum[2*N]; 
  newnum[0] = ntip+1; //root is after tips
  int numtip = 1;     //tips are first
  int numint = ntip+2; //internal nodes are after the root
  for( j=1;j<2*N;j++){//extant species first
    if(edge_alive[j-1]){
      newnum[j]=numtip;
      numtip++;
    }
  }
  for( j=1;j<2*N;j++){
    if(!edge_alive[j-1]){
      if( is_tip[j] ){//extinct species are second
	newnum[j]=numtip;
	numtip++;
      }else{//internal nodes are last
	newnum[j]=numint;
	numint++;
      }
    }
  }
  for( j=0;j<is*4;j++){
    P_edge[j] = newnum[P_edge[j]];
  }

  for( j=0;j<NH-1;j++){
    free(nlistt[j]);
    free(tcur[j]);
  }
  free(nlistt);
  free(tcur);
  return ntip;
}


int main(int argc, char ** argv)//takes the path to an input file as argument
{
  /*----------------Parameters--------------------*/
  //records the time of execution
  time_t start,end;         //beginning and ending times
  float dif;                //difference between beginning and ending

  //model parameters
  // default values
  int numbsim = 1;
  int nh = 219; //observed number of hosts
  int Nmin = 50;
  float lp = 0.96;  //P spe rate
  float mp = 0.59;  //P extinction rate
  float pTmrca = 13.96;
  char* fin   = "";
  char* fout  = "";
  char* fileP = "distrib.txt";
   
  float pco = -1.0;
  float hs  = -1.0;

  time (&start);//records the beginning time
  int ss = time(NULL);
  int c;

  // input
  opterr = 0;
  while ((c = getopt (argc, argv, "l:m:c:s:t:n:N:P:o:i:S:?")) != -1)
    switch (c)
      {
      case 'l':
        lp = atof(optarg);
        break;
      case 'm':
        mp = atof(optarg);
        break;
      case 'c':
        pco = atof(optarg);
        break;
      case 's':
        hs = atof(optarg);
        break;
      case 't':
        pTmrca = atof(optarg);
        break;
      case 'n':
        numbsim = atoi(optarg);
        break;
      case 'N':
	Nmin = atoi(optarg);
        break;
      case 'P':
        fileP = optarg;
        break;
      case 'o':
        fout = optarg;
        break;
      case 'i':
        fin = optarg;
        break;
      case 'S':
        ss = atof(optarg);
        break;
      case '?':
        if (optopt == 'e')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }
  
  srand(ss); 
  
  //parameters
  float pcol[numbsim]; //cospeciation probabilities
  float lswitchl[numbsim]; //host switch rates
    
  //read host phylogeny
  int i,j;
  FILE * fr;
  FILE * fr2;
  char foedge[100];
  char foedgelength[100];
  sprintf(foedge, "%shost_edge.txt",fin);
  sprintf(foedgelength, "%shost_edgelength.txt",fin);
  
  printf("Input files: %s, %s\n",foedge,foedgelength);
  fr  = fopen(foedge,"rb");   
  int nout,nout2;
  nout = fscanf(fr,"%d",&(nh));
  nh --;
  printf("nh=%d\n",nh);
  fclose(fr);

  int nhe = 2*(nh-1); //observed number of edges 
  int   host_E[nhe*2]; //edge list
  float host_E_length[nhe]; //edge length

  fr  = fopen(foedge,"rb");   
  fr2 = fopen(foedgelength,"rb");   
  for( i=0; i<nhe; i++){
    nout = fscanf(fr,"%d %d",&(host_E[2*i]),&(host_E[2*i+1]));
    nout2 = fscanf(fr2,"%f",&(host_E_length[i]));
  }
  fclose(fr);
  fclose(fr2);
  
  printf("Input host phylogeny:\n");
  for( i=0; i<nhe; i++){
    nout = printf("edge %d: %d -> %d,\t",i,(host_E[2*i]),(host_E[2*i+1]));
    nout2 = printf("length=%f\n",(host_E_length[i]));
  }
  
  int nfP = 1;
  nout = 0;
  if(fileP!=NULL){
    FILE * fP;
    fP  = fopen(fileP,"rb");   
    nout = 1;
    while ( (c=fgetc(fP)) != EOF ) {
      if( (c==' ')|(c=='\t') ) nfP++;
    }
    fclose(fP);
  }else{
    nfP = 61;
  }
  float dge[nfP]; //distribution of number of H per P (dge[0]: specialists, dge[i>0]: generalists)

  if(fileP!=NULL){
    int ntmp;
    FILE * fP;
    fP  = fopen(fileP,"rb");   
    printf("distribution of Parasite per Host: \t");
    for( i=0; i<nfP; i++){
      fscanf(fP,"%f",&(dge[i]));
      printf("%f ",dge[i]);
    }
    printf("\n");
    fclose(fP);
  }else{//default value
    dge[0]=0.456621; dge[1]=0.193660527776115; dge[2]=0.08881349485094; dge[3]=0.0463596318170004; dge[4]=0.0285872373542903; dge[5]=0.020633127585547; dge[6]=0.0166342494882219; dge[7]=0.0142732101971074;   dge[8]=0.0126298786420637; dge[9]=0.0113351384170491; dge[10]=0.0102375180693405; dge[11]=0.00927183136964164; dge[12]=0.00840738741099738;  dge[13]=0.00762754523217738; dge[14]=0.00692161829188089;  dge[15]=0.00628164706554527; dge[16]=0.00570109261294602;  dge[17]=0.00517429000670136; dge[18]=0.00469620400075332;  dge[19]=0.00426230642496276; dge[20]=0.00386850395980884;  dge[21]=0.00351108796542828; dge[22]=0.00318669500417025;  dge[23]=0.00289227339660586; dge[24]=0.0026250538072553;   dge[25]=0.00238252291807734; dge[26]=0.00216239968007697;  dge[27]=0.00196261381646778; dge[28]=0.00178128633456101;  dge[29]=0.00161671184708468; dge[30]=0.0014673425303588;   dge[31]=0.00133177356596966; dge[32]=0.00120872992804381;  dge[33]=0.00109705439147095; dge[34]=0.000995696648141144; dge[35]=0.000903703428775051;dge[36]=0.00082020953741926;  dge[37]=0.00074442971427707; dge[38]=0.000675651250340947; dge[39]=0.000613227284365852;dge[40]=0.00055657071914096;  dge[41]=0.000505148699842076;dge[42]=0.000458477602533591; dge[43]=0.00041611848568682; dge[44]=0.000377672961936251; dge[45]=0.000342779451247578;dge[46]=0.000311109780258572; dge[47]=0.000282366095809603;dge[48]=0.000256278063635581; dge[49]=0.000232600325873009;dge[50]=0.000211110193470021; dge[51]=0.000191605551796526;dge[52]=0.000173902959756724; dge[53]=0.000157835924526156;dge[54]=0.000143253335687192; dge[55]=0.000130018044036018;dge[56]=0.000118005571694784; dge[57]=0.000107102941397543;dge[58]=0.0000972076139394068; dge[59]=0.0000882265237956343;dge[60]=0.0000800752038406542; //from fitted distribution of mixture model
  }
  float Pg = dge[0];//proportion of specialist P 
 

  /***************  main loop ******************/
  int nmax = 10000; //maximum size of phylogeny above which simulation is stopped; ensures that we do not get unrealistically large phylogenies which make computation times explode

  printf("Parasite speciation rate:=%f, Parasite extinction rate=%f, Parasite TMRCA=%f\n",lp,mp,pTmrca);
    
  for( i=0;i<numbsim;i++){//numbsim simulations
    //create result tables
    int   res_E[nmax*2*2]; //edge list
    for( j=0;j<nmax*2*2;j++) res_E[j] = 0;
    float res_E_length[nmax*2]; //edge length
    for( j=0;j<nmax*2;j++) res_E_length[j] = 0;

    int** pedges = malloc( nmax*2*sizeof(int*)); // parasite position (edge)
    for( j=0;j<nmax*2;j++) pedges[j] = malloc(nh*sizeof(int));
    int nedges[nmax*2]; // parasite position (edge)
    for( j=0;j<nmax*2;j++) nedges[j] = 0;
    bool edge_alive[nmax*2];
    for( j=0;j<nmax*2;j++) edge_alive[j] = false;
    //random parameters
    if(pco<0) pcol[i]     = rand()/(float)RAND_MAX; 
    else pcol[i] = pco;
    if(hs<0) lswitchl[i] = rand()/(float)RAND_MAX; 
    else lswitchl[i] = hs;
    printf("- loop %i: proba cospe=%f, host-switch rate=%f\n",i,pcol[i],lswitchl[i]);
    //run simulation
    int N = bd_P_cospe(host_E,host_E_length,pedges,nedges,edge_alive,nh,lp,mp,1-Pg,nfP,dge,lswitchl[i],pcol[i],pTmrca,res_E,res_E_length,nmax);
    int nalive=0;
    for( j=0;j<2*N;j++){
      if(edge_alive[j]) nalive++;
    }

    int attempts = 0;
    const int max_attempts = 2000; // avoids infinite loops when parameters are too "extinction-heavy"
    while(nalive < Nmin && attempts < max_attempts){
      attempts++;
      N = bd_P_cospe(host_E,host_E_length,pedges,nedges,edge_alive,nh,lp,mp,1-Pg,nfP,dge,lswitchl[i],pcol[i],pTmrca,res_E,res_E_length,nmax);
      printf("\tNumber of tips: %d\n", N);

      nalive = 0;
      for( j=0;j<2*N;j++){
        if(edge_alive[j]) nalive++;
      }
    }

    if(nalive < Nmin){
      fprintf(stderr, "WARNING: Could not reach Nmin=%d alive tips after %d attempts (last N=%d, nalive=%d). Consider increasing -t (Tmrca), increasing -l, decreasing -m, or lowering -N.\n",
              Nmin, max_attempts, N, nalive);
      // Continue anyway with the last simulated tree (may be empty). If you prefer to abort, replace the next line with: return 1;
    }
    printf("\tend compute cophylogeny\n");
    //print results
    char fedge[100];
    char fhpassoc[100];
    char fedgelength[100];
    char fnnode[100];
    sprintf(fedge, "%sedges_%d.txt",fout,i);
    sprintf(fhpassoc, "%shpassoc_%d.txt",fout,i);
    sprintf(fedgelength, "%sedgelength_%d.txt",fout,i);
    sprintf(fnnode, "%snnode_%d.txt",fout,i);
    FILE * fe;
    FILE * ff;
    FILE * fg;
    FILE * fh;
    fe=fopen(fedge,"wb");  
    ff=fopen(fhpassoc,"wb");  
    fg=fopen(fedgelength,"wb");  
    fh=fopen(fnnode,"wb");  
    for( j=0;j<(N-1)*2;j++){
      fprintf(fe,"%d\t%d\n",res_E[j*2], res_E[j*2+1]);
      fprintf(fg,"%.30f\n",res_E_length[j]);
    }
    fprintf(fh,"%d\n",N-1);
    printf("\tend write edges\n");
    int icur = 1;
    int iii = 0;
    while (icur <= nalive) {
      if (edge_alive[iii]) {
        if (nedges[iii + 1] == 0) printf("ERROR nedges\n");
        for (j = 0; j < nedges[iii + 1]; j++) {
          int hedge = pedges[iii + 1][j];
          if (hedge < 0 || hedge >= nhe) {
            fprintf(stderr, "WARNING: invalid host edge index %d for parasite_tip=%d (iii=%d, j=%d). Skipping association.\n",
                    hedge, icur, iii, j);
            continue;
          }
          fprintf(ff, "%d\t%d\n", icur + 2 * nh - 1, host_E[2 * hedge + 1]);
        }
        icur++;
      }
      iii++;
    }
    printf("\tend write associations\n");
    fclose(fe);
    fclose(ff);
    fclose(fg);
    fclose(fh);
    
    //write list of parameters (drawn from uniform prior distribution)
    char fparam[100];
    sprintf(fparam, "%sparams_%d.txt",fout,i);
    FILE * fi;
    fi=fopen(fparam,"wb");  
    fprintf(fi,"%lf\t%lf\n",pcol[i], lswitchl[i]);
    fclose(fi);

    //delte tables
    for( j=0;j<nmax*2;j++){
      free(pedges[j]);
    }
    free(pedges);
  }

  
  time (&end);//records ending time
  dif = difftime (end,start);
  printf ("Finished. It took  %.2lf min\n", dif/60.0 );

}
