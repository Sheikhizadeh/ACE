#include <string>
#include <sys/types.h>
#include <unistd.h> 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <bitset>

#ifdef __linux__
#include <malloc.h>
#include <sys/sysinfo.h>
#endif

#ifdef __APPLE__
#include <sys/sysctl.h>
#define malloc_trim(x)	{}
#endif


using namespace std;

class Read
{
public:
  short length;
  unsigned short bytes[MAXREADLEN/8+1];
};


//  The class for nodes of the K_mer trie. Each node has 4 children and one parent. 
class node
{
  public: node* children[4];
  node* parent;
  node(node* p=0)
  {
    parent=p;
    for(int i=0;i<4;++i)
      children[i]=0;
  }
  ~node()
  {
    for(int i=0;i<4;++i)
      if(children[i])
      {
	delete children[i];
	children[i]=0;
      }
  }
};

// Some global variables
#ifdef __linux__
struct sysinfo info;
#endif

double teta,Maxteta,cov,decrs;
char sym[]={ 'A', 'C', 'G', 'T' };
int p1_len,p2_len,p_len,p1_num,p_num,phase,K, k, cores, top,phnum,minL=100000,maxL=0;
long avgL;
long size;
unsigned long N, G;
Read* Reads;
node*** root;
node** level;
node** prev;
bool* first;
unsigned long Index[20];
char filetype;

// To convert the integer prefix to the character_string prefix 
char* code(char* pc, int prefix)
{
  int i;
  for(i=0;i<p1_len;++i)
    pc[i]=sym[(prefix>>2*(p1_len-i-1)) & 3];
  return pc;
}

// To extract a long from a char array starting with some whitespaces 
long atol(char* line)
{
  int i,j;
  for(i=0;!isdigit(line[i]);++i);
  for(j=0;isdigit(line[i]);++i,++j)
    line[j]=line[i];
  line[j]=0;
  return strtoul(line,0,10);
}

// To execute a shell command and returning the results 
void exec(char* cmd,char* result) 
{
  FILE* pipe = popen(cmd, "r");
  if (!pipe) 
  {
    cerr<<"Unable to execute the shell command."<<endl<<flush;
    exit(1);
  }
  char buffer[128];
  result[0]=0;
  while(!feof(pipe)) 
    if(fgets(buffer, 128, pipe) != NULL)
      strcat(result,buffer);
    pclose(pipe);
}

// To insert a k_mer compressed in the form of a long integer into the trie.
// "r" is the number of the read to whom the k_mer belongs, and "loc" its location in that read.
// "x" is the number of the subtir which the k_mer should be inserted in.
void insert(unsigned long kmer, unsigned long r,  short loc, short x)
{
  int i , myK=K,t=top,p=p_len,j,number;
  node *v;
  kmer=kmer>>(2*p);
  for (number = 0, i = p+1; i <= p+t; ++i,kmer=kmer>>2 )
    number = number * 4 + (kmer & 3L);
  if(root[x][number]==0)
    root[x][number]=new node();
  v = root[x][number];
  for (; i <= myK; ++i,kmer=kmer>>2 )
  {
    j = (kmer & 3L);
    if (!v->children[j])
      v->children[j] = new node(v);
    v = v->children[j];
  }
  if (v->children[0] == 0)
  {
    if (first[x])
    {
      first[x] = false;
      prev[x] = level[x] = v;
    }
    else
    {
      prev[x]->children[3] = v;
      prev[x] = v;
      v->children[3] = 0;
    }
  }
  v->children[0]=(node *)(((long)v->children[0])+1);
  v->children[1]=(node *)(r);
  v->children[2]=(node *)((long)loc);
}

// To create the subtrie which will contain all K_mers beginning with "prefix".  
// "prefix" is the integer equivalent of the real character string prefix, e.g. 0 for AA, 1 for AC and 15 for TT. 
// The subtries of this subtires are constructed in parallel, by "cores" omp threads.
void create_trie(int prefix)
{
  for (int i = 0; i < cores; ++i)
    first[i] = true;
  omp_set_num_threads(cores);
  #pragma omp parallel
  {
    short loc,j,myK=K,y;
    unsigned long tmp,fp,rp,number,frame=(1L<<(2*p_len))-1,i, myN=N, begin,set=~(3L<<2*K);
    short L,p=p_len;
    y=omp_get_thread_num();
    begin=N/cores*y+1;
    for(j=0,number=0,tmp=y;j<p2_len;++j,tmp=tmp>>2)
      number=number*4+(tmp&3L);
    for(j=0,    tmp=prefix;j<p1_len;++j,tmp=tmp>>2)
      number=number*4+(tmp&3L);
    for (i = begin; i <= myN; ++i)
    {
      L=Reads[i].length;
      for (j = 0, fp=rp = 0; j<myK ; ++j)
      {
	fp =(fp<<2) |   ((Reads[i].bytes[(myK-j-1) / 8] >> (2 * (7-(myK-j-1) % 8)) ) & 3L);
	rp =(rp<<2) |(3-((Reads[i].bytes[(      j) / 8] >> (2 * (7-(j      ) % 8)) ) & 3L));
      }
      for (loc = 0; loc < L-myK+1; ++loc)
      {
	if((fp & frame) == number)
	  insert(fp, i, loc, y);
	if((rp & frame) == number)
	  insert(rp, i, -(loc+myK-1), y);
	fp= (fp>>2)|(((Reads[i].bytes[(loc+myK) / 8] >> (2 * (7-(loc+myK) % 8))) & 3L) <<(2*(myK-1)));
	rp=((rp<<2) & set ) | (3-((Reads[i].bytes[(loc+myK) / 8] >> (2 * (7-(loc+myK) % 8))) & 3L));
      }
    }
    for (i = 1; i < begin; ++i)
    {
      L=Reads[i].length;
      for (j = 0, fp=rp = 0; j<myK ; ++j)
      {
	fp =(fp<<2) |   ((Reads[i].bytes[(myK-j-1) / 8] >> (2 * (7-(myK-j-1) % 8)) ) & 3L);
	rp =(rp<<2) |(3-((Reads[i].bytes[(      j) / 8] >> (2 * (7-(j      ) % 8)) ) & 3L));
      }
      for (loc = 0; loc < L-myK+1; ++loc)
      {
	if((fp & frame) == number)
	  insert(fp, i, loc, y);
	if((rp & frame) == number)
	  insert(rp, i, -(loc+myK-1), y);
	fp= (fp>>2)|(((Reads[i].bytes[(loc+myK) / 8] >> (2 * (7-(loc+myK) % 8))) & 3L) <<(2*(myK-1)));
	rp=((rp<<2) & set ) | (3-((Reads[i].bytes[(loc+myK) / 8] >> (2 * (7-(loc+myK) % 8))) & 3L));
      }
    }
  }
}

// To detect and correct the substitutions in a parallel manner. "cores" threads scan the leaves' chain seeking for low_frequent k_mers 
void detect_and_correct(int pfx)
{
  
  omp_set_num_threads(cores);
  #pragma omp parallel
  {
    int y=omp_get_thread_num();
    bool found;
    node* hpoint;
    node* branch;
    node* scan;
    node* next;
    int i, j, k, d, l1, l2, l3,myK=K,min=max(top,(K-p_len)/2)+p_len+1;
    unsigned short x=3;
    double foot = teta, maxfoot = Maxteta;
    long r, loc,location;
    int s[100];
    for (scan = level[y]; scan != 0; scan = scan->children[3])
    {
      if ( ((long)scan->children[0]) > 0 && ((long)scan->children[0]) <= foot)
      {
	found = false;
	branch = scan;
	for (i = myK; i >= min && !found; --i)
	{
	  for(j=0;j<4;++j)
	    if (branch->parent->children[j] == branch)
	    {
	      s[i] = j;
	      break;
	    }
	    branch = branch->parent;
	  for (l1 = 1,s[i] = (s[i] + 1) % 4; l1 <= 3 && !found; ++l1,s[i] = (s[i] + 1) % 4)
	  {
	    hpoint = branch;
	    for (d = i; d <= myK && hpoint; ++d)
	      hpoint = hpoint->children[ s[d] ];
	    if (hpoint && (((long)hpoint->children[0]) > maxfoot ) )
	    {
	      found = true;
	      location=(long)scan->children[2];
	      r = (long)scan->children[1];
	      if(location>=0)
	      {
		loc=location + i-1 ;
		#pragma omp critical (c1) 
		{
		  Reads[r].bytes[loc / 8] &= ( ~(x<<(2 * (7-loc % 8))));
		  Reads[r].bytes[loc / 8] |= (s[i]<<(2 * (7-loc % 8)));
		}
	      }
	      else
	      {
		loc=-location - i+1 ;
		#pragma omp critical (c2) 
		{
		  Reads[r].bytes[loc / 8] &= (~(x<<(2 * (7-loc % 8))));
		  Reads[r].bytes[loc / 8] |= ((3-s[i])<<(2 * (7-loc % 8)));
		}
	      }
	    }//if fits
	  }//for l
	}//for i
	//  2-change search just lunches from second round on because its very time-consuming to run it on the first round where the bank contains many errors                        
	if(phase>1) 
	{                        
	  
	  branch = scan;
	  for (i = myK ; i >= min; --i)
	  {
	    for(j=0;j<4;++j)
	      if (branch->parent->children[j] == branch)
	      {
		s[i] = j;
		break;
	      }
	      branch = branch->parent;
	  }
	  for (i = myK; i >= min+1 && !found; --i)
	    for (l1 = 1,s[i] = (s[i] + 1) % 4; l1 <= 3 && !found; ++l1,s[i] = (s[i] + 1) % 4)
	    {
	      for(j=i-1;j>=min && !found; --j)
		for (l2 = 1,s[j] = (s[j] + 1) % 4; l2 <= 3 && !found; ++l2,s[j] = (s[j] + 1) % 4)
		{				
		  hpoint = branch;
		  for (d = min; d <= myK && hpoint; ++d)
		    hpoint = hpoint->children[ s[d] ];
		  if (hpoint && ((long)(hpoint->children[0]) > maxfoot ) )
		  {
		    found = true;
		    location=(long)scan->children[2];
		    r = (long)scan->children[1];
		    if(location>=0)
		    {
		      loc=location + i-1;
		      #pragma omp critical (c3) 
		      {
			Reads[r].bytes[loc / 8] &= ( ~(x<<(2 * (7-loc % 8))));
			Reads[r].bytes[loc / 8] |= (s[i]<<(2 * (7-loc % 8)));
		      }
		      loc=location + j-1;
		      #pragma omp critical (c4) 
		      {
			Reads[r].bytes[loc / 8] &= ( ~(x<<(2 * (7-loc % 8))));
			Reads[r].bytes[loc / 8] |= (s[j]<<(2 * (7-loc % 8)));
		      }
		    }
		    else
		    {
		      loc=-location - i+1;
		      #pragma omp critical (c5) 
		      {
			Reads[r].bytes[loc / 8] &= (~(x<<(2 * (7-loc % 8))));
			Reads[r].bytes[loc / 8] |= ((3-s[i])<<(2 * (7-loc % 8)));
		      }
		      loc=-location - j+1;
		      #pragma omp critical (c6) 
		      {			  
			Reads[r].bytes[loc / 8] &= (~(x<<(2 * (7-loc % 8))));
			Reads[r].bytes[loc / 8] |= ((3-s[j])<<(2 * (7-loc % 8)));
		      }
		    }
		  }//if fits
		}//for l2
	    }//for l 
	}//if phase
      }// if scan
    }//for scan
    for(scan=level[y];scan!=0;scan=next)
    {
      next=scan->children[3];
      scan->children[0]=0;
      scan->children[1]=0;
      scan->children[2]=0;
      scan->children[3]=0;
    }	
  }
}

// To free the memory occupied by the subtires of the K_mer trie
void free_trie()
{
  int i,j;
  for (i = 0; i < cores; ++i)
    for (j = 0; j < size; ++j)
      if(root[i][j])
      {
	delete root[i][j];
	root[i][j]=0;
      }
}

// To write the corrected reads into the output file 
void write_reads(char* argv[])
{
  int i,j,f, k, l;
  unsigned short x=3;
  char ch;
  char s[10000];
  cout<<"Storing the corrected archive(s)...\n"<<flush;
  ifstream input(argv[2]);
  if(!input )
  {
    cout<<"Failed to open "<<argv[2]<<flush;
    exit(1);
  }
  ofstream output(argv[3]);
  if(!output )
  {
    cout<<"Failed to open "<<argv[3]<<flush;
    exit(1);
  }
  if(filetype=='@')
    for (i=1; i <= N; ++i)
    {
      input.getline(s,10000,'\n');
      output<<s<<endl;
      input.getline(s,10000,'\n');
      for (j = 0; j<Reads[i].length; ++j)
	output<<sym[ (Reads[i].bytes[j / 8] >> (2 * (7-j % 8) )) & x ];
      output<<endl;
      input.getline(s,10000,'\n');
      output<<s<<endl;
      input.getline(s,10000,'\n');
      output<<s<<endl;
    }
    else
      for (i=1,input.getline(s,10000,'\n'); i <= N; ++i)
      {
	output<<s<<endl; 
	for(k=0,input.getline(s,10000,'\n');s[0]!=filetype && ! input.eof();input.getline(s,10000,'\n'),output<<endl)
	  for (j = 0; s[j]!=0 && s[j]!='\r'; ++j,++k)
	    output<<sym[(Reads[i].bytes[k / 8] >> (2 * (7-k % 8) )) & x];;
      }
      input.close();
    output.close();
}

// To delete all the allocated memory and return it to OS before finishing the run
void finalize(char* argv[])
{
  int i,j;
  write_reads(argv);
  cout<<"Finalizing the run ...\n"<<flush;
  delete Reads;
  for (i = 0; i < cores; ++i)
  {
    for (j = 0; j < size; ++j)
      if(root[i][j])
      {
	delete root[i][j];
	root[i][j]=0;
      }
  }	    
  delete root;
  delete level;
  delete prev;
  delete first;
}

// To load reads from the input file and store them in a binary compressed format
void load_reads(char* filename)
{
  char command[10000];
  char s[10000];
  char line[1000];
  char out[10000];
  int f,j,k,l,c,d;
  unsigned long i;
  unsigned short x;
  cout<<"Loading reads ...\r"<<flush;
  ifstream input(filename);
  if(!input )
  {
    cout<<"Failed to open "<<filename<<endl<<flush;
    exit(1);
  }
  input.get(filetype);
  input.close();
  if(filetype=='@')
  {
    strcpy(command,"wc -l ");
    strcat(command,filename);
    exec(command,out);
    N += atol(out)/4;
  }
  else
  {
    strcpy(command,"grep '>' ");
    strcat(strcat(command,filename),"| wc -l");
    exec(command, out);
    N += atol(out);
  }
  Reads=new Read[N+1];
  input.open(filename);
  input.getline(s,10000,'\n');
  if(filetype=='@')
    for(i=1;!input.eof();++i)
    {
      input.getline(s,10000,'\n');
      l=strlen(s);
      avgL+=l;
      if(l<minL)
	minL=l;
      if(maxL<l)
	maxL=l;
      d=(int)ceil(l/8);
      Reads[i].length=l;
      for(j=0;j<=d;++j)
	Reads[i].bytes[j]=0;
      for (j = 0; s[j]!=0 && s[j]!='\r'; ++j)
      {
	x = ((s[j]=='N'?rand()%4:Index[s[j]-65]) << (2 * (7 - j % 8)));
	Reads[i].bytes[j / 8] |= x;
      }
      if(i%1000000==0)
	cout<<"Loading reads ..."<<"("<<i/1000000<<"M)\r"<<flush;
      input.getline(s,10000,'\n');
      input.getline(s,10000,'\n');
      input.getline(s,10000,'\n');
    }
    else
      for(i=1;!input.eof();++i)
      {
	for(s[0]=0,input.getline(line,100,'\n');line[0]!='>' && ! input.eof();input.getline(line,100,'\n'))
	  strcat(s,line);
	l=strlen(s);
	avgL+=l;	  
	if(l<minL)
	  minL=l;
	if(maxL<l)
	  maxL=l;
	d=(int)ceil(l/8);
	Reads[i].length=l;
	for(j=0;j<=d;++j)
	  Reads[i].bytes[j]=0;
	for (j = 0; s[j]!=0 && s[j]!='\r'; ++j)
	{
	  x = ((s[j]=='N'?rand()%4:Index[s[j]-65]) << (2 * (7 - j % 8)));
	  Reads[i].bytes[j / 8] |= x;
	}
	if(i%1000000==0)
	  cout<<"Loading reads ..."<<"("<<i/1000000<<"M)\r"<<flush;	  
      }
      input.close();
      avgL/=N;
}

// To initialize the environment
void Initialize(char* argv[])
{
  int i,j;
  long n;
  long long freemem;
  G = strtoul(argv[1],0,10);
  N=0;
  Index[0] = 0;
  Index[2] = 1;
  Index[6] = 2;
  Index[19] = 3;
  p2_len=(int)min(2.0,floor(log2(sysconf(_SC_NPROCESSORS_ONLN))/2));
  cores=1<<2*p2_len;
  root=new node**[cores];
  level=new node*[cores];
  prev=new node*[cores];
  first=new bool[cores];	
  
  #ifdef __linux__
  sysinfo(&info);
  freemem = 1+info.totalram/1000000000*info.mem_unit;
  #endif
  
  #ifdef __APPLE__
  int32_t freepages, specpages;
  size_t len = sizeof(freepages);
  sysctlbyname("vm.page_free_count", &freepages, &len, NULL, 0);
  sysctlbyname("vm.page_speculative_count", &specpages, &len, NULL, 0);
  freemem = 1+(freepages+specpages)/1000000000*4096;
  #endif
  
  p1_len=max(2.0,3.0-log2(freemem)/2);
  p1_num=1<<2*p1_len;
  p_len=p1_len+p2_len;
  p_num=1<<2*p_len;
  phnum = 7;
  cout<<"=====================   ACE   ====================\n"<<flush;
  load_reads(argv[2]);
  K = floor(log2(G)/2)+10;
  k = K - p1_len;
  cov = N * (float)avgL  / G;
  top =(int)(min(floor(log2(N)/2),floor(log2(G)/2))- p_len);
  size = 1L<<(2*top);
  teta = 2+cov/20;
  Maxteta = 2*teta;
  decrs = phnum>1?(teta-2) / (phnum - 1):0;
  for (i = 0; i < cores; ++i)
  {
    root[i] = new node*[size];
    for (j = 0; j < size; ++j)
      root[i][j]=0;
  }
  cout<<"Genome length = "<<G<<endl;
  cout<<"Reads length  = "<<minL<<" to "<<maxL<<endl;
  cout<<"Reads count   = "<<N<<endl;
  cout<<"Coverage      = "<<cov<<endl;
  cout<<"K             = "<<K<<endl;
  cout<<"............................\n"<<flush;
}

//,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, The main() ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
int main(int argc, char* argv[])
{
  time_t ptime=0,ctime;
  time(&ptime); 
  if(argc<4)
  {
    cout<<"--------------------------------------------------------------------\n";
    cout<<"ACE corrects substitution errors in a single/paired_end Illumina archive.\n";
    cout<<"It receives the genome length and the name of input file and output files.\n";
    cout<<"To run ACE:\n";
    cout<<"./ace Genome_length(bp) Inputfile Outputfile\n\n"<<flush;
    return 0;
  }
  Initialize(argv);
  time(&ctime);
  cout<<(ctime-ptime)<<" seconds elapsed.              \n\n"<<flush;
  int i;
  int pfx;
  char pc[10];
  pc[p1_len]=0;
  for (phase = 1; phase <= phnum ; ++phase)
  {
    cout<< "Round "<<phase<<" of "<<phnum<<" : teta = "<<(int)teta<<endl<<flush;
    for (pfx = 0; pfx < p1_num ; ++pfx)
    {  
      cout<<"Subtrie "<<code(pc,pfx)<<" : constructing... "<<flush;
      create_trie(pfx);
      cout<<"processing... "<<flush;
      detect_and_correct(pfx);
      cout<<"done.\n"<<flush;

      #ifdef __linux__
      sysinfo(&info);
      if(info.freeram < 2*info.totalram/3) 
      #endif

      #ifdef __APPLE__
      char buf[100];
      size_t buflen = 100;
      long hwmem, freepages, specpages;
      sysctlbyname("hw.memsize", &buf, &buflen, NULL, 0);
      hwmem = atol(buf);
      sysctlbyname("vm.page_free_count", &buf, &buflen, NULL, 0);
      freepages = atol(buf);
      sysctlbyname("vm.page_speculative_count", &buf, &buflen, NULL, 0);
      specpages = atol(buf);
      if ((freepages+specpages)*4096 < hwmem/2) 
      #endif

      {
	free_trie();
	malloc_trim(0);
      }
    }//pfx
    time(&ctime);
    cout<<(ctime-ptime)<<" seconds elapsed.              \n\n"<<flush;
    teta    -= decrs;
    Maxteta =2*teta;
  }//phase
  finalize(argv);
  time(&ctime);
  cout<< "\nTotal Runtime : "<<ctime-ptime<< " seconds\n";
  char line[100];
  sprintf(line, "//proc//%d//status", getpid());
  ifstream file(line);
  for(file.getline(line,100,'\n');strncmp(line,"VmPeak:",7)!=0;file.getline(line,100,'\n'));
  cout<< "Maximum space : "<<atol(line)/1024<<" Mb"<<endl;
  cout<<"===================== Success ====================\n"<<flush;
  return 1;
}
