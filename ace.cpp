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
// The class for storing each short read in memory, which contains:
// seq: The character string of the read
// length: Length of the read in characters
// locations: an array in which location of K_mers beggining with a specific prefix is stored, i.e. 
// those who start with AA stored in indexes 0 to ends[0], AC in indexes ends[0]+1 to ends[1], and ...
class Read
{
  public: char* seq;
  public: unsigned char length;
  public: short int* locations;
  public: short int* ends;
  public: Read()
  {
    seq=0;
    length=0;
    ends=0;
    locations=0;
  }
  public: Read(int L,int minK, int p2_num)
  {
    seq=new char[L+1];
    seq[0]=0;
    length=L;
    ends = new short int[p2_num];
    locations = new  short int[2*(length - minK + 1)];
  }
  public: ~Read()
  {
    delete seq;
    delete ends;
    delete locations;
  }
};
// Some global variables
#ifdef __linux__
struct sysinfo info;
#endif

int Ks[8];
double Foot,MaxFoot,cov,decrs,incrs;
char sym[]={ 'A', 'C', 'G', 'T' };
int p1_len,p2_len,p1_num,p2_num,phase,K, k, minK, cr, cores, L,top,phnum;
double total=0,count=0;
long int size;
long int N=0, G;
Read **Reads;
Read **reads;
node*** root;
node** level;
node** prev;
bool* first;
unsigned long Index[20];
short int*** temp;
int** indx;
char fasta;
char *infile1=0, *infile2=0, *outfile1=0, *outfile2=0; 
// This function makes the arrays "location" and "ends" of a read, when they have been loaded or have been changed.
// "int r" is the number of read and "int c" is the number of subtires which should be constructed in parallel. 
void make(int r, int c)
{
  
  int i, j, number;
  short int loc, k = -1;
  short int**tmp=temp[c];
  int* inx=indx[c];
  for(j=0;j<p2_num;++j)
    inx[j]=0;
  for (loc = 1; loc <= L-minK+1; ++loc)
  {
    for (j = 0, number = 0; j < p2_len; ++j)
      number =number*4+Reads[r]->seq[loc+j];
    tmp[number][inx[number]] = loc;
    inx[number]++;
  }
  for (loc = L; loc >= minK; --loc)
  {
    for (j = 0, number = 0; j < p2_len; ++j)
      number = number * 4 + (3-Reads[r]->seq[loc - j]);
    tmp[number][inx[number]] = -loc;
    inx[number]++;
  }
  for (i = 0; i < p2_num; ++i)
  {
    for (j = 0; j < inx[i]; ++j)
      Reads[r]->locations[++k] = tmp[i][j];
    Reads[r]->ends[i] = k;
  }
}
// This function inserts the forward spelling of a k_mer into the trie.
// "r" is the number of the read to whom the k_mer belongs, and "loc" its location in that read.
// "x" is the number of the subtir which the k_mer should be inserted in.
void insertF(long r,  long loc, int x)
{
  int i = 0, j, myk=k,t=top,myp=p1_len;
  int number;
  node *v;
  for (number = 0, i = cr+1; i < t; ++i)
    number = number * 4 + Reads[r]->seq[loc + myp + i - 1];
  //number %= size;
  if(root[x][number]==0)
    root[x][number]=new node();
  v = root[x][number];
  for (; i <= myk; ++i)
  {
    j = Reads[r]->seq[loc + myp + i - 1];
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
  v->children[2]=(node *)(loc);
}
// This function inserts the reverse-complement spelling of a k_mer into the trie.
// "r" is the number of the read to whom the k_mer belongs, and "loc" its location in that read.
// "x" is the number of the subtir which the k_mer should be inserted in. 
void insertR(long r,  long loc, int x)
{
  int i = 0, j, myk=k,t=top,myp=p1_len;
  int number;
  node *v;
  for (number = 0, i = cr+1; i < t; ++i)
    number = number * 4 + (3 - Reads[r]->seq[-loc - myp - i + 1]);
  //number %= size;
  if(root[x][number]==0)
    root[x][number]=new node();
  v = root[x][number];
  for (; i <= myk; ++i)
  {
    j = (3 - Reads[r]->seq[-loc - myp - i + 1]) ;
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
  v->children[2]=(node *)(loc);
}
// This function creates the subtrie which will contain all K_mers beginning with "prefix".  
// "prefix" is the integer equivalent of the real character string prefix, e.g. 0 for AA, 1 for AC and 15 for TT. 
// The subtries of this subtires are constructed in parallel, by "cores" omp threads.
void create_trie(int prefix)
{
  for (int i = 0; i < cores; ++i)
    first[i] = true;
  omp_set_num_threads(cores);
  #pragma omp parallel
  {
    int l,begin,end, y;
    long loc,i;	    
    y=omp_get_thread_num();
    for (i = 1; i <= N; ++i)
    {
      begin = (y == 0 && prefix == 0 ? 0 : Reads[i]->ends[cores * prefix + y - 1] + 1);
      end = Reads[i]->ends[cores * prefix + y];
      for (l = begin; l <= end; ++l)
      {
	loc = Reads[i]->locations[l];
	if ((loc > 0 && loc+K-1<=L))
	  insertF(i, loc, y);
	else if((loc < 0 && -loc-K+1>=1 ) )
	  insertR(i, loc, y);
      }
    }
  }
}
// This function detects and corrects the substitutions in a parallel manner. "cores" threads scan the leaves' chain
// seeking for low_frequent k_mers,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 
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
    int i, j, d, l1, l2, tmp1,tmp2,myk=k,myp=p1_len,t=top;
    double foot = Foot, maxfoot = MaxFoot;
    long r, loc,location;
    int* s = new int[myk];
    for (scan = level[y]; scan != 0; scan = scan->children[3])
    {
      if ( ((long)scan->children[0]) > 0 && ((long)scan->children[0]) <= foot)
      {
	found = false;
	branch = scan;
	for (i = myk - 1; i >= t && !found; --i)
	{
	  for(j=0;j<4;++j)
	    if (branch->parent->children[j] == branch)
	    {
	      s[i] = j;
	      break;
	    }
	    branch = branch->parent;
	  for (l1 = 1; l1 <= 3 && !found; ++l1)
	  {
	    tmp1 = s[i];
	    s[i] = (s[i] + l1) % 4;
	    hpoint = branch;
	    for (d = i; d < myk && hpoint; ++d)
	      hpoint = hpoint->children[ s[d] ];
	    if (hpoint && (((long)hpoint->children[0]) > maxfoot ) )
	    {
	      found = true;
	      location=(long)scan->children[2];
	      r = (long)scan->children[1];
	      if(location>0)
	      {
		loc=location + myp + i;
		if(Reads[r]->seq[loc]!=s[i])
		{
		  #pragma omp critical (c1) 
		  {
		    Reads[r]->seq[loc] = s[i];
		    Reads[r]->seq[0] = phase ;
		  }
		}
	      }
	      else
	      {
		loc=-location - myp - i;
		if(Reads[r]->seq[loc]!=3-s[i])
		{
		  #pragma omp critical (c2) 
		  {
		    Reads[r]->seq[loc] = 3-s[i];
		    Reads[r]->seq[0] = phase ;
		  }
		}
	      }
	    }//if fits
	    s[i] = tmp1;
	  }//for l
	}//for i
	//  2-change search just lunches from second round on 
	//  because its very time-consuming to run it on the first round where the bank contains many errors                        
	if(phase>1) 
	{                        
	  
	  branch = scan;
	  for (i = myk - 1; i >= top; --i)
	  {
	    for(j=0;j<4;++j)
	      if (branch->parent->children[j] == branch)
	      {
		s[i] = j;
		break;
	      }
	      branch = branch->parent;
	  }
	  for (i = myk - 1; i >= top && !found; --i)
	    for (l1 = 1; l1 <= 3 && !found; ++l1)
	    {
	      tmp1 = s[i];
	      s[i] = (s[i] + l1) % 4;
	      for(j=i-1;j>=top && !found; --j)
		for (l2 = 1; l2 <= 3 && !found; ++l2)
		{				
		  tmp2 = s[j];
		  s[j] = (s[j] + l2) % 4;
		  hpoint = branch;
		  for (d = top; d < myk && hpoint; ++d)
		    hpoint = hpoint->children[ s[d] ];
		  if (hpoint && ((long)(hpoint->children[0]) > maxfoot ) )
		  {
		    found = true;
		    location=(long)scan->children[2];
		    r = (long)scan->children[1];
		    if(location>0)
		    {
		      loc=location + myp + i;
		      if(Reads[r]->seq[loc]!=s[i])
		      {
			#pragma omp critical (c3) 
			{
			  Reads[r]->seq[loc] = s[i];
			  Reads[r]->seq[0] = phase ;
			}
		      }
		      loc=location + myp + j;
		      if(Reads[r]->seq[loc]!=s[j])
		      {
			#pragma omp critical (c4) 
			{
			  Reads[r]->seq[loc] = s[j];
			  Reads[r]->seq[0] = phase ;
			}
		      }
		    }
		    else
		    {
		      loc=-location - myp - i;
		      if(Reads[r]->seq[loc]!=3-s[i])
		      {
			#pragma omp critical (c5) 
			{
			  Reads[r]->seq[loc] = 3-s[i];
			  Reads[r]->seq[0] = phase ;
			}
		      }
		      loc=-location - myp - j;
		      if(Reads[r]->seq[loc]!=3-s[j])
		      {
			#pragma omp critical (c6) 
			{
			  Reads[r]->seq[loc] = 3-s[j];
			  Reads[r]->seq[0] = phase ;
			}
		      }
		    }
		  }//if fits
		  s[j]=tmp2;
		}//for l2
		s[i] = tmp1;
	    }//for l
	}
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
  #pragma omp parallel
  {
    int y=omp_get_thread_num(),begin,end,r;
    begin= y*N/cores+1 ;
    end= y==cores-1 ? N :(y+1)*N/cores ;
    for(r=begin;r<=end;++r)
      if(Reads[r]->seq[0] == phase)
      {
	#pragma omp critical (c7) 
	{
	  ++count;
	}	
	make(r,y);
      }
  }    
  
}
// ,,,,,,,,,,,,,,,,,,,This function frees the memory occupied by the subtires of the K_mer trie,,,,,,,,,,
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
// ,,,,,,,this function deleted all the allocated memory to return it to OS before finishing the run,,,,,,
void finalize()
{
  int i,j;
  for(i=1;i<=N;++i)
    delete Reads[i];
  for (i = 0; i < cores; ++i)
  {
    for (j = 0; j < size; ++j)
      if(root[i][j])
      {
	delete root[i][j];
	root[i][j]=0;
      }
      for(j=0;j<p2_num;++j)
	delete temp[i][j];
      delete temp[i];
    delete indx[i];
  }	    
  delete temp;
  delete indx;
  delete root;
  delete level;
  delete prev;
  delete first;
}
// ,,,,,,,,,,,,,,,,,,,,,,,,,This function writes the corrected reads into new file(s),,,,,,,,,,,,,,,,,,,,, 
void correctBank(int i, int n, char* infile, char* outfile)
{
  int j, k, l;
  char ch;
  char s[10000];
  ifstream input(infile);
  if(!input )
  {
    cout<<"Failed to open "<<infile<<flush;
    exit(0);
  }
  ofstream output(outfile);
  if(!output )
  {
    cout<<"Failed to open "<<outfile<<flush;
    exit(0);
  }
  if(fasta=='@')
    for (input.getline(s,10000,'\n'); i <= n; ++i)
    {
      output<<s<<endl; 
      for(l=k=0,input.getline(s,10000,'\n');s[0]!='+' && ! input.eof();input.getline(s,10000,'\n'),++l,output<<endl)
	for (j = 0; s[j]!=0 && s[j]!='\r'; ++j,++k)
	  output<<sym[Reads[i]->seq[k+1]];
	output<<s<<endl; 
      for(input.getline(s,10000,'\n');l>0;input.getline(s,10000,'\n'),--l)
	output<<s<<endl; 
    }
    else
      for (input.getline(s,10000,'\n'); i <= n; ++i)
      {
	output<<s<<endl; 
	for(k=0,input.getline(s,10000,'\n');s[0]!=fasta && ! input.eof();input.getline(s,10000,'\n'),output<<endl)
	  for (j = 0; s[j]!=0 && s[j]!='\r'; ++j,++k)
	    output<<sym[Reads[i]->seq[k+1]];
      }
      input.close();
    output.close();
}
//,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,To obtain the length of reads in filename ,,,,,,,,,,,,,,,,,,,,,,,,,,,,
int reads_length(char* filename)
{
  char s[10000];
  int i,j,k,l=0,n=N;
  ifstream input(filename);
  if(!input )
  {
    cout<<"Failed to open "<<filename<<endl<<flush;
    exit(1);
  }
  input.getline(s,10000,'\n');
  fasta=s[0];
  if(fasta=='@')
    for(i=0,input.getline(s,10000,'\n');s[0]!='+' ;input.getline(s,10000,'\n'),++i) 
      for (j = 0; s[j]!=0 && s[j]!='\r'; ++j)
	l++;
      else
	for(i=0,input.getline(s,10000,'\n');s[0]!='>' ;input.getline(s,10000,'\n'),++i) 
	  for (j = 0; s[j]!=0 && s[j]!='\r'; ++j)
	    l++;
	  input.close();
	return l;
}
//,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,To load reads of filename into the memory ,,,,,,,,,,,,,,,,,,,,,,,,,,,,
int load_reads(char* filename)
{
  char s[10000];
  int i,j,k,l;
  if(filename!=0)
  {
    ifstream input(filename);
    if(!input )
    {
      cout<<"Failed to open "<<filename<<endl<<flush;
      exit(1);
    }	  
    input.getline(s,10000,'\n');
    fasta=s[0];
    if(fasta=='@')
      for(i=1;!input.eof();++i)
      {
	Reads[N+i]=new Read(L, minK, p2_num);
	for(l=k=0,input.getline(s,10000,'\n');s[0]!='+' && ! input.eof();input.getline(s,10000,'\n'),++l)
	  for (j = 0; s[j]!=0 && s[j]!='\r'; ++j,++k)
	    Reads[N+i]->seq[k+1]=(s[j]=='N'?rand()%4:Index[s[j]-65]);
	  for(input.getline(s,10000,'\n');l>0;input.getline(s,10000,'\n'),--l);
	  make(N+i,0);
      }
      else
	for(i=1;!input.eof();++i)
	{
	  Reads[N+i]=new Read(L, minK, p2_num);
	  for(k=0,input.getline(s,10000,'\n');s[0]!=fasta && ! input.eof();input.getline(s,10000,'\n'))
	    for (j = 0; s[j]!=0 && s[j]!='\r'; ++j,++k)
	      Reads[N+i]->seq[k+1]=(s[j]=='N'?rand()%4:Index[s[j]-65]);
	    make(N+i,0);
	}
	N+=i-1;
      input.close();
  }
  return i-1;
}
// ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,To initialize the environment,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
void Initialize()
{
  int i,j;
  long n;
  float freemem;
  N=0;
  Index[0] = 0;
  Index[2] = 1;
  Index[6] = 2;
  Index[19] = 3;
  cr=(int)floor(log2(sysconf(_SC_NPROCESSORS_ONLN))/2);
  cores=1<<2*cr;
  root=new node**[cores];
  level=new node*[cores];
  prev=new node*[cores];
  first=new bool[cores];	

  #ifdef __linux__
  sysinfo(&info);
  freemem = (float) info.freeram*(float) info.mem_unit/1000000000.0;
  #endif

  #ifdef __APPLE__
  int32_t freepages, specpages;
  size_t len = sizeof(freepages);
  sysctlbyname("vm.page_free_count", &freepages, &len, NULL, 0);
  sysctlbyname("vm.page_speculative_count", &specpages, &len, NULL, 0);
  freemem = (((float) freepages+(float) specpages)*4096.0)/1000000000.0;
  #endif
  
  p1_len=max(2.0,4-log2(freemem)/2);
  p1_num=1<<2*p1_len;
  p2_len=p1_len+cr;
  p2_num=1<<2*p2_len;
  phnum = 7;
  Reads=new Read*[500000000];
  temp=new short int**[cores];
  indx=new int*[cores];
  L=reads_length(infile1);
  for (i = 0; i < cores; ++i)
  {	    
    temp[i]=new short int*[p2_num];
    indx[i]=new int[p2_num];
    for(j=0;j<p2_num;++j)
    {
      temp[i][j]= new short int[ 2 * (L - minK + 1)];
      indx[i][j]=0;
    }
  }
  cout<<"Loading reads...\n"<<flush;
  cout<<load_reads(infile1)<<" reads loaded from "<<infile1<<endl<<flush;
  if(infile2)
    cout<<load_reads(infile2)<<" reads loaded from "<<infile2<<endl<<flush;
  K = floor(log2(N)/2)+floor(log2(G)/2);
  Ks[1]=K;
  Ks[2]=K+2;
  Ks[3]=K+1;
  Ks[4]=K;
  Ks[5]=K-2;
  Ks[6]=K-1;
  Ks[7]=K;
  minK=K-2;
  Read** pr=Reads;
  Reads=new Read*[N+1];
  for(i=1;i<=N;++i)
    Reads[i]=pr[i];
  delete pr;
  malloc_trim(0);
  cov = N * (float)L  / G;
  top = (int)floor(log2(G)/2)-1;
  size = 1<<(2*(top-1));
  Foot = 2 + cov / 10;
  MaxFoot = 3 + cov / 10;
  decrs = phnum>1?(Foot-2) / (phnum - 1):0;
  for (i = 0; i < cores; ++i)
  {
    root[i] = new node*[size];
    for (j = 0; j < size; ++j)
      root[i][j]=0;
  }
  cout<<"Genome = " << G << "bp\n"<< "Read = " << L << "bp\nCoverage = " << cov << "\nPrefix length = "<<p1_len<< "\n...............................\n\n"<<flush;
}
//,,,,,,,,,,,,,,,,,,,,,,,,,,,,, To convert the integer prefix to the character_string prefix ,,,,,,,,,,,,,,,,
char* code(char* pc, int prefix)
{
  int i;
  for(i=0;i<p1_len;++i)
    pc[i]=sym[(prefix>>2*(p1_len-i-1)) & 3];
  return pc;
}
//,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, To add "corrected" to the name of input file .........,,,
char* newname(char* name)
{
  int l= strlen(name);
  char c[]="corrected";
  char* s=new char[l+11];
  int i,j;
  strcpy(s,name);
  for(i=l; s[i]!='.'; --i)
    s[i+10]=s[i];
  s[i+10]='.';
  for(i+=9,j=8;j>=0;--i,--j)
    s[i]=c[j];
  return s;
}

//,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, The main() ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
int main(int argc, char* argv[])
{
  int i;
  if(argc>2)
  {
    G = strtoull(argv[1],0,10);
    infile1 = argv[2];
    outfile1 = newname(infile1);
    if(argc==4)
    {
      infile2 = argv[3];
      outfile2 = newname(infile2);
    }
  }
  else
  {
    cout<<"--------------------------------------------------------------------\n";
    cout<<"ACE corrects substitution errors in single/paired Illumina archives.\n";
    cout<<"To run ACE:\n";
    cout<<"./ace Genome_length(bp) Inputfile(s)\n\n"<<flush;
    return 0;
  }
  Initialize();
  time_t ptime=0,ctime;
  time(&ptime); 
  int pfx;
  char* pc=new char[p1_len+1];
  pc[p1_len]=0;
  for (phase = 1; phase <= phnum ; ++phase)
  {
    K=Ks[phase];
    k = K - p1_len;
    cout<< "Round "<<phase<<" of "<<phnum<<":\tK="<< K<<", teta="<<floor(Foot)<<endl<<flush;
    for (count=0,pfx = 0; pfx < p1_num ; ++pfx)
    {  
      cout<<"Processing "<<code(pc,pfx)<<" subtrie ...\r"<<flush;
      create_trie(pfx);
      detect_and_correct(pfx);
      #ifdef __linux__
      sysinfo(&info);
      if(info.freeram < info.totalram/5)
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
      if ((freepages+specpages)*4096 < hwmem/5)
      #endif
      {
	free_trie();
	malloc_trim(0);
      }
    }//pfx
    free_trie();
    malloc_trim(0);
    time(&ctime);
    total+=count;
    cout<<ctime-ptime<<" seconds elapsed.\t\n"<<fixed<<count<<" corrections to reach "<<fixed<<total<<" corrections.\n\n"<<flush;
    Foot   -= decrs;
    MaxFoot-= decrs/2;
  }//phase
  time(&ctime);
  cout<< "Total runtime : "<<ctime-ptime<< " s.\n";
  cout<<total<<" bases corrected.\n"<<flush;
  if(infile2)
  {
    correctBank(1,N/2,infile1,outfile1);
    correctBank(N/2+1,N,infile2,outfile2);
  }
  else
    correctBank(1,N,infile1,outfile1);
  finalize();
  return 1;
}
