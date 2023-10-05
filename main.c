//Bootstrap RMSE
//by D. Ashley Robinson


//Includes
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>


//Defines
#define             kVERSION                    "version 1.0"
#define             kDATE                       "September 20, 2022"
#define             kMAXSITES                   1380 //should be 1 longer than needed as seqs stored as null-terminated strings
#define             kMAXPOPULATIONS             2
#define             kMAXISOLATESPERPOPULATION   154
#define             kMAXISOLATES                (kMAXPOPULATIONS*kMAXISOLATESPERPOPULATION)


//Globals
FILE                *infile, *outfile;

unsigned int        seed;
int                 numSamples, numBoots, numPermutations;
int                 BootCounter;
int                 numIsolates, numSites, numPopulations;
int                 numIsolatesPerPopulation[kMAXPOPULATIONS], samplenumIsolatesPerPopulation[kMAXPOPULATIONS];
int                 firstIsolatePerPopulation[kMAXPOPULATIONS], samplefirstIsolatePerPopulation[kMAXPOPULATIONS];//array position of 1st isolate from theIsolates
int                 derivedAlleleCounter[kMAXPOPULATIONS][kMAXSITES];
double              derivedAlleleFreqs[kMAXPOPULATIONS][kMAXSITES];
double              SSE;//sum of squared errors

char                theIsolatesName[kMAXISOLATES][20];//need to implement check on taxon labels >20 chars, these labels not used atm
char                theIsolatesData[kMAXISOLATES][kMAXSITES], sampletheIsolatesData[kMAXISOLATES][kMAXSITES];


//Prototypes
void DoReadInfile(void);
void DoRunIt(void);
void DoAnalyze(void);
void DoCleanupSample(void);
void DoBootSample(void);
void DoPermutation(void);
void DoOpenTheFile(void);
void DoCloseTheFile(void);


//Main
int main(void){
    int     procFlg;
    char    inname[50];
    
    printf("-----------------------------------------------------\n");
    printf("Bootstrap RMSE\n");
    printf("%s, %s\n", kVERSION, kDATE);
    printf("by D. Ashley Robinson (darobinson@umc.edu)\n");
    printf("Department of Cell and Molecular Biology\n");
    printf("University of Mississippi Medical Center\n");
    printf("Jackson, MS 39216\n");
    printf("-----------------------------------------------------\n\n");
    
    printf("Enter name of infile: ");
    scanf("%s", inname);
    infile=fopen(inname,"r");
    DoReadInfile();
    fclose(infile);
    
    procFlg=0;
    do{
        printf("\nEnter pseudorandom number seed [Range=1-65534]: ");
        scanf("%d", &seed);
        if((seed>=1)&&(seed<=65534)){procFlg=1;}
    }while(procFlg==0);
    srand(seed);
    
    printf("\nEnter per population sample size to draw: ");
    scanf("%d", &numSamples);

    printf("\nEnter number of bootstrap replicates: ");
    scanf("%d", &numBoots);
    
    printf("\nCalculating...\n");
    
    DoRunIt();
    
    printf("\nCalculations complete.\n");
    
    return(0);
}


//DoReadInfile
void DoReadInfile(void){
    int     i, j, counter;
    
    if(infile==NULL){
        printf("\nThere is an error with the infile.\n");
        exit(0);
    }
    
    //read parameters
    fscanf(infile, "%d", &numIsolates);
    if(numIsolates>kMAXISOLATES){
        printf("\nToo many isolates (maximum=%d).\n", kMAXISOLATES);
        exit(0);
    }
    
    fscanf(infile, "%d", &numSites);
    if(numSites>kMAXSITES){
        printf("\nToo many sites (maximum=%d).\n", kMAXSITES);
        exit(0);
    }
    
    fscanf(infile, "%d", &numPopulations);
    if(numPopulations>kMAXPOPULATIONS){
        printf("\nToo many populations (maximum=%d).\n", kMAXPOPULATIONS);
        exit(0);
    }
    
    counter=0;
    for(i=0; i<numPopulations; i++){
        fscanf(infile, "%d", &numIsolatesPerPopulation[i]);
        samplenumIsolatesPerPopulation[i]=numIsolatesPerPopulation[i];
        counter+=numIsolatesPerPopulation[i];
        if(numIsolatesPerPopulation[i]>kMAXISOLATESPERPOPULATION){
            printf("\nToo many isolates in population %d (maximum=%d).\n", i, kMAXISOLATESPERPOPULATION);
            exit(0);
        }
    }
    if(numIsolates!=counter){
        printf("\nTotal number of isolates (%d) does not match sum of isolates in subpopulations (%d).\n", numIsolates, counter);
    }
    
    //determine index of first isolate in each population, 0 to n-1 range with n as numisolates
    counter=1;
    for(i=0; i<numPopulations; i++){
        firstIsolatePerPopulation[i]+=counter-1;
        samplefirstIsolatePerPopulation[i]=firstIsolatePerPopulation[i];
        counter+=numIsolatesPerPopulation[i];
    }
    
    //read taxonlabels and sequences
    //reading both as strings
    //can printf both as strings eg %s %s taxa[i] seqs[i]
    //or reference specific sites in seq array and print as char eg %c seqs[i][j]
    
    for(i=0; i<numIsolates; i++){
        fscanf(infile, "%s", theIsolatesName[i]);
        fscanf(infile, "%s", theIsolatesData[i]);
    }
    
    for(i=0; i<numIsolates; i++){
        for(j=0; j<numSites; j++){
            sampletheIsolatesData[i][j]=theIsolatesData[i][j];
        }
    }
    
}


//DoRunIt
//master function
void DoRunIt(void){
    int     i;
    
    DoOpenTheFile();
    
    BootCounter=0;
    DoAnalyze();
    
    for(i=0; i<numBoots; i++){
        BootCounter=i;
        DoCleanupSample();
        DoBootSample();
        DoAnalyze();
    }
        
    DoCloseTheFile();
    
}


//DoAnalyze
void DoAnalyze(void){
    int s, i, j, k, l;
    double c1, n1, d1;
    
    //first make allele counts for each site and population
    for(s=0; s<numSites; s++){
        
        for(i=0; i<numPopulations; i++){
            j=samplefirstIsolatePerPopulation[i];
            k=(samplefirstIsolatePerPopulation[i]+samplenumIsolatesPerPopulation[i]);
            
            for(l=j; l<k; l++){
                if(sampletheIsolatesData[l][s]=='1'){derivedAlleleCounter[i][s]++;}
            }
        }
    }
    
    //check
    //for(i=0; i<numSites; i++){
    //    for(j=0; j<numPopulations; j++){
    //       fprintf(outfile, "%d ", derivedAlleleCounter[j][i]);
    //    }
    //    fprintf(outfile, "\n");
    //}

    //calcfreqs
    for(i=0; i<numSites; i++){
        for(j=0; j<numPopulations; j++){
            c1=(double)derivedAlleleCounter[j][i];
            n1=(double)samplenumIsolatesPerPopulation[j];
            derivedAlleleFreqs[j][i]=(c1/n1);
            //fprintf(outfile, "%.5f ", derivedAlleleFreqs[j][i]);
        }
        
        d1=(derivedAlleleFreqs[0][i]-derivedAlleleFreqs[1][i]);
        SSE+=(d1*d1);
        //check
        //fprintf(outfile, "%.9f ", SSE);
        //fprintf(outfile, "\n");
    }
    fprintf(outfile, "Bootcount %d\t%.9f\n", BootCounter, (SSE/numSites));
    
}


//DoCleanupSample
void DoCleanupSample(void){
    int     i, j;
    
    SSE=0;
    for(i=0; i<numPopulations; i++){
        for(j=0; j<numSites; j++){
            derivedAlleleCounter[i][j]=0;
            derivedAlleleFreqs[i][j]=0;
        }
    }
    
    
}


//DoBootSample
void DoBootSample(void){
    int r, s, m, n, i, j, procFlg;
    
    m=0;
    for(i=0; i<numPopulations; i++){
        
        n=0;
        procFlg=0;
        samplenumIsolatesPerPopulation[i]=numSamples;
        samplefirstIsolatePerPopulation[i]=m;
        
        do{
            r=((unsigned int) rand()%numIsolatesPerPopulation[i]);//eg rand()%10 gives 0-9
            s=(firstIsolatePerPopulation[i]+r);
            //check
            //fprintf(outfile, "%d ", s);
            for(j=0; j<numSites; j++){
                sampletheIsolatesData[m][j]=theIsolatesData[s][j];
            }
            m++;
            n++;
            if(n==numSamples){procFlg=1;}
            
        }while(procFlg==0);
        
        
        //fprintf(outfile, "\n");
    }
    
    //check
    //for(i=0; i<numIsolates; i++){
    //    fprintf(outfile, "%s", sampletheIsolatesData[i]);
    //   fprintf(outfile, "\n");
    //}
    
}


//DoOpenTheFile
void DoOpenTheFile(void){
    char    name[50];
    
    sprintf(name, "outfile.txt");
    outfile=fopen(name, "w");
    strcpy(name, "");
}


//DoCloseTheFile
void DoCloseTheFile(void){
    fclose(outfile);
}

