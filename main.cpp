//  Individual base simulation coding the basic model of mite mating strategy evolution.
//  Haploid model, no assexual reproduction; With male precedence: females give new offspring only for their first mating event. Males and females return to mating pool after sex. Male guards forever.
//  We define two genes: x for guarding and y as a neutral gene.
//  May 2017. Code with VARYING POPULATION SIZE, where mortality decreases population size and reproduction replenishes it.

//PARTIAL DISCRIMINATION
//Change the way matings are decided.


#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <cassert>
using namespace std;


gsl_rng *rng_r; // gnu scientific rng


const int Npop = 10000; //Population size
double muFjf, muFjp, muFsf, muFsm, muFmf, muFmm, muMgf, muMgp, muMsf, muMsm, mortall, Tf, Tm, Tg, g, mf, mm, cg, cs, Mg, Ms, sexratio1, minim, maxim, mu, stdev,k, discr, assex, prop;
int init_uFjf,init_uFjp, init_uFsf, init_uFsm, init_uFmf, init_uFmm, init_uMgf, init_uMgp, init_uMsf, init_uMsm, indivF, indivM, Nsumpop;
int NuFjf,NuFjp, NuFsf, NuFsm, NuFmf, NuFmm, NuMgf, NuMgp, NuMsf, NuMsm; // counters for all four classes
int tmax, skip, points, repmax, sum_events, switcher, events;
double Tfo, Tmo;
unsigned long seed = 123454321; //seed for random number generation
double xcount, ycount, meanx, meany, stdevx, stdevy; //to calculate mean and stdev
unsigned int nprobs; //large number for the total of possible events
double init_x;
double init_y;
double init_xs;
double init_ys;
bool do_stats=1;
double prob;

#define memsize 200


//  Individual struct: used to distinguish each individual within its class.
struct Individual
{
    double x;   //  guarding gene
    double y;   //  neutral gene
    double xs;  //  sperm guarding gene
    double ys;  //  sperm neutral gene
};

typedef Individual Population[Npop];
Population FEMAjuvfree;
Population FEMAjuvpaired;
Population FEMAsexfree;
Population FEMAsexmate;
Population FEMAmaturefree;
Population FEMAmaturemate;
Population MALEguardfree;
Population MALEguardpaired;
Population MALEsearchfree;
Population MALEsexmate;

//Initial arguments (eventually to be uploaded in invocation)
void initArguments(int argc, char *argv[])
{
    
    //    //Input parameters to run in the shell
        int j=0;
        repmax = atoi(argv[++j]);    //number of repetitons
        tmax = atoi(argv[++j]); //max time per repetition
        mu = atof(argv[++j]);  //mutation rate
        stdev = atof(argv[++j]);   //standard deviation of mutational effects
        mortall = atof(argv[++j]); //mortality
        cs = atof(argv[++j]); //cost of searching
        cg = atof(argv[++j]); //cost of guarding
        Ms = atof(argv[++j]);    //effectiveness of search encounters
        Mg = atof(argv[++j]);    //effectiveness of guarding encounters
        Tf = atoi(argv[++j]);   // Time to sexual maturation of juvenile females
        Tg = atoi(argv[++j]);   //Time guarding
        Tfo = atoi(argv[++j]);  //Time-out of females after mating   
        Tmo = atoi(argv[++j]);  //Time-out of males after mating
        sexratio1 = atof(argv[++j]);     //primary sex ratio
        discr = atof(argv[++j]);    //searcher male descrimination rate
        assex = atof(argv[++j]); // probability of assexual reproduction
        seed = atoi(argv[++j]); //seed
        prop = atof(argv[++j]); //initial proportions of male classes
    
        if (argc!=19)
        {
            cout  <<  "start the program like this:\n" << argv[0]
            << " <repmax> <tmax> <mut rate> <stdev mut> <mortality>  <cs> <cg> <Ms> <Mg> <Tf> <Tg> <Tfo> <Tmo> <sexratio1> <discrim rate> <assex rate> <seed> <prop> \n"
            << endl;
            exit (-1);
        }
        cout<<"Input parameters: "<<repmax<<" "<<tmax<<" "<<" "<<mu<<" "<<" "<<stdev<<" "<<" "<<mortall<<" "<<cs<<" "<<cg<<" "<<" "<<Ms<<" "<<Mg<<" "<<Tf<<" "<<Tg<<" "<<Tfo<<" "<<Tmo<<" "<<sexratio1<<" "<<" "<<discr<<" "<<assex<<" "<<" "<<prop<<" "<<seed<<endl;
    
    
    muFjf = muFjp = muFsf = muFsm = muFmf = muFmm = muMgf = muMgp = muMsf = muMsm = mortall;    
    
    //Initial numbers for each class.
    init_uFjf = Npop/2;
    init_uFjp = 0;
    init_uFsf = 0;
    init_uFsm = 0;
    init_uFmf = 0;
    init_uFmm = 0;
    init_uMgf = (Npop/2)*prop;
    init_uMgp = 0;
    init_uMsf = (Npop/2)*(1-prop);
    init_uMsm = 0;
    
    //Initial chromosome values
    init_x = 0.5;
    init_y = 0.5;
    
    //To determine when we print values
    skip = 100000;
}

//Function to initialize class values
void Init()
{
    
    //Initial class values:
    NuFjf = init_uFjf;
    NuFjp = init_uFjp;
    NuFsf = init_uFsf;
    NuFsm = init_uFsm;
    NuFmf = init_uFmf;
    NuFmm = init_uFmm;
    NuMgf = init_uMgf;
    NuMgp = init_uMgp;
    NuMsf = init_uMsf;
    NuMsm = init_uMsm;
    
    
    
    //Initial trait values within each class
    for (int i=0 ; i < NuFjf ; i++)
    {
        FEMAjuvfree[i].x = init_x;
        FEMAjuvfree[i].y = init_y;
        FEMAjuvfree[i].xs = init_x;
        FEMAjuvfree[i].ys = init_y;
    }
    
    for (int i=0 ; i < NuFjp ; i++)
    {
        FEMAjuvpaired[i].x = init_x;
        FEMAjuvpaired[i].y = init_y;
        FEMAjuvpaired[i].xs = init_x;
        FEMAjuvpaired[i].ys = init_y;
    }
    
    for (int i=0 ; i < NuFsf ; i++)
    {
        FEMAsexfree[i].x = init_x;
        FEMAsexfree[i].y = init_y;
        FEMAsexfree[i].xs = init_x;
        FEMAsexfree[i].ys = init_y;
    }
    
    for (int i=0 ; i < NuFsm ; i++)
    {
        FEMAsexmate[i].x = init_x;
        FEMAsexmate[i].y = init_y;
        FEMAsexmate[i].xs = init_x;
        FEMAsexmate[i].ys = init_y;
    }
    
    for (int i=0 ; i < NuFmf ; i++)
    {
        FEMAmaturefree[i].x = init_x;
        FEMAmaturefree[i].y = init_y;
        FEMAmaturefree[i].xs = init_x;
        FEMAmaturefree[i].ys = init_y;
    }
    
    for (int i=0 ; i < NuFmm ; i++)
    {
        FEMAmaturemate[i].x = init_x;
        FEMAmaturemate[i].y = init_y;
        FEMAmaturemate[i].xs = init_x;
        FEMAmaturemate[i].ys = init_y;
    }
    
    for (int i=0 ; i < NuMgf ; i++)
    {
        MALEguardfree[i].x = init_x;
        MALEguardfree[i].y = init_y;
        MALEguardfree[i].xs = init_x;
        MALEguardfree[i].ys = init_y;
    }
    
    for (int i=0 ; i < NuMgp ; i++)
    {
        MALEguardpaired[i].x = init_x;
        MALEguardpaired[i].y = init_y;
        MALEguardpaired[i].xs = init_x;
        MALEguardpaired[i].ys = init_y;
    }
    
    for (int i=0 ; i < NuMsf ; i++)
    {
        MALEsearchfree[i].x = init_x;
        MALEsearchfree[i].y = init_y;
        MALEsearchfree[i].xs = init_x;
        MALEsearchfree[i].ys = init_y;
    }
    
    for (int i=0 ; i < NuMsm ; i++)
    {
        MALEsexmate[i].x = init_x;
        MALEsexmate[i].y = init_y;
        MALEsexmate[i].xs = init_x;
        MALEsexmate[i].ys = init_y;
    }
}


//Mortality functions for each class
int death_FEMAjuvfree(Individual *FEMAjuvfree, int indiv_index, int clasnumbs)
{
    //    cout<<" death_FEMAjuvfree" <<endl;
    assert(clasnumbs > 0);
    FEMAjuvfree[indiv_index].x = FEMAjuvfree[NuFjf - 1].x;
    FEMAjuvfree[indiv_index].y = FEMAjuvfree[NuFjf - 1].y;
    FEMAjuvfree[indiv_index].xs = FEMAjuvfree[NuFjf - 1].xs;
    FEMAjuvfree[indiv_index].ys = FEMAjuvfree[NuFjf - 1].ys;
    assert(FEMAjuvfree[indiv_index].x > 0 && FEMAjuvfree[indiv_index].y > 0);
    
    clasnumbs--;
    return clasnumbs;
}

int death_FEMAjuvpaired(Individual *death_FEMAjuvpaired, int indiv_index, int clasnumbs)
{
    //    cout<<" death_FEMAjuvpaired" <<endl;
    assert(clasnumbs > 0);
    death_FEMAjuvpaired[indiv_index].x = death_FEMAjuvpaired[NuFjp - 1].x;
    death_FEMAjuvpaired[indiv_index].y = death_FEMAjuvpaired[NuFjp - 1].y;
    death_FEMAjuvpaired[indiv_index].xs = death_FEMAjuvpaired[NuFjp - 1].xs;
    death_FEMAjuvpaired[indiv_index].ys = death_FEMAjuvpaired[NuFjp - 1].ys;
    clasnumbs--;
    return clasnumbs;
}

int death_FEMAsexfree(Individual *FEMAsexfree, int indiv_index, int clasnumbs)
{
    //    cout<<" death_FEMAsexfree" <<endl;
    assert(clasnumbs > 0);
    FEMAsexfree[indiv_index].x = FEMAsexfree[NuFsf - 1].x;
    FEMAsexfree[indiv_index].y = FEMAsexfree[NuFsf - 1].y;
    FEMAsexfree[indiv_index].xs = FEMAsexfree[NuFsf - 1].xs;
    FEMAsexfree[indiv_index].ys = FEMAsexfree[NuFsf - 1].ys;
    clasnumbs--;
    return clasnumbs;
}

int death_FEMAsexmate(Individual *FEMAsexmate, int indiv_index, int clasnumbs)
{
    //    cout<<" death_FEMAsexmate" <<endl;
    assert(clasnumbs > 0);
    FEMAsexmate[indiv_index].x = FEMAsexmate[NuFsm - 1].x;
    FEMAsexmate[indiv_index].y = FEMAsexmate[NuFsm - 1].y;
    FEMAsexmate[indiv_index].xs = FEMAsexmate[NuFsm - 1].xs;
    FEMAsexmate[indiv_index].ys = FEMAsexmate[NuFsm - 1].ys;
    clasnumbs--;
    return clasnumbs;
}

int death_FEMAmaturefree(Individual *FEMAmaturefree, int indiv_index, int clasnumbs)
{
    //    cout<<" death_FEMAmaturefree" <<endl;
    assert(clasnumbs > 0);
    FEMAmaturefree[indiv_index].x = FEMAmaturefree[NuFmf - 1].x;
    FEMAmaturefree[indiv_index].y = FEMAmaturefree[NuFmf - 1].y;
    FEMAmaturefree[indiv_index].xs = FEMAmaturefree[NuFmf - 1].xs;
    FEMAmaturefree[indiv_index].ys = FEMAmaturefree[NuFmf - 1].ys;
    clasnumbs--;
    return clasnumbs;
}
int death_FEMAmaturemate(Individual *FEMAmaturemate, int indiv_index, int clasnumbs)
{
    //    cout<<" death_FEMAmaturemate" <<endl;
    assert(clasnumbs > 0);
    FEMAmaturemate[indiv_index].x = FEMAmaturemate[NuFmm - 1].x;
    FEMAmaturemate[indiv_index].y = FEMAmaturemate[NuFmm - 1].y;
    FEMAmaturemate[indiv_index].xs = FEMAmaturemate[NuFmm - 1].xs;
    FEMAmaturemate[indiv_index].ys = FEMAmaturemate[NuFmm - 1].ys;
    clasnumbs--;
    return clasnumbs;
}

int death_MALEguardfree(Individual *MALEguardfree, int indiv_index, int clasnumbs)
{
    //    cout<<" death_MALEguardfree" <<endl;
    assert(clasnumbs > 0);
    MALEguardfree[indiv_index].x = MALEguardfree[NuMgf - 1].x;
    MALEguardfree[indiv_index].y = MALEguardfree[NuMgf - 1].y;
    MALEguardfree[indiv_index].xs = MALEguardfree[NuMgf - 1].xs;
    MALEguardfree[indiv_index].ys = MALEguardfree[NuMgf - 1].ys;
    clasnumbs--;
    return clasnumbs;
}

int death_MALEguardpaired(Individual * MALEguardpaired, int indiv_index, int clasnumbs)
{
    //    cout<<" death_MALEguardpaired" <<endl;
    assert(clasnumbs > 0);
    MALEguardpaired[indiv_index].x = MALEguardpaired[NuMgp - 1].x;
    MALEguardpaired[indiv_index].y = MALEguardpaired[NuMgp - 1].y;
    MALEguardpaired[indiv_index].xs = MALEguardpaired[NuMgp - 1].xs;
    MALEguardpaired[indiv_index].ys = MALEguardpaired[NuMgp - 1].ys;
    clasnumbs--;
    return clasnumbs;
}

int death_MALEsearchfree(Individual *MALEsearchfree, int indiv_index, int clasnumbs)
{
    //    cout<<" death_MALEsearchfree" <<endl;
    assert(clasnumbs > 0);
    MALEsearchfree[indiv_index].x = MALEsearchfree[NuMsf - 1].x;
    MALEsearchfree[indiv_index].y = MALEsearchfree[NuMsf - 1].y;
    MALEsearchfree[indiv_index].xs = MALEsearchfree[NuMsf - 1].xs;
    MALEsearchfree[indiv_index].ys = MALEsearchfree[NuMsf - 1].ys;
    clasnumbs--;
    return clasnumbs;
}

int death_MALEsexmate(Individual *MALEsexmate, int indiv_index, int clasnumbs)
{
    //    cout<<" death_MALEsexmate" <<endl;
    assert(clasnumbs > 0);
    MALEsexmate[indiv_index].x = MALEsexmate[NuMsm - 1].x;
    MALEsexmate[indiv_index].y = MALEsexmate[NuMsm - 1].y;
    MALEsexmate[indiv_index].xs = MALEsexmate[NuMsm - 1].xs;
    MALEsexmate[indiv_index].ys = MALEsexmate[NuMsm - 1].ys;
    clasnumbs--;
    return clasnumbs;
}


//Mutation function
//  We limit genes to have a value between 0 and 1 (this is different from Kokko's model) - see             lines 372, 473, 906
void mutate(double &mX, double &mY)
{
    mX += gsl_rng_uniform(rng_r)<mu ? gsl_ran_gaussian(rng_r, stdev) : 0;
    mY += gsl_rng_uniform(rng_r)<mu ? gsl_ran_gaussian(rng_r, stdev) : 0;
    
    mX = mX < 0 ? 0.000001 : mX;
    mX = mX > 1 ? 1 : mX;
    
    mY = mY < 0 ? 0.000001 : mY;
    mY = mY > 1 ? 1 : mY;
    
}

//Breeding of new individual function
void create_kid(Individual &mother, Individual &kid)
{
    assert(mother.x > 0 && mother.y > 0 && mother.xs > 0 && mother.ys > 0);
    
    kid.x = gsl_rng_uniform(rng_r) < 0.5 ? mother.x : mother.xs;
    kid.y = gsl_rng_uniform(rng_r) < 0.5 ? mother.y : mother.ys;
    
    assert(kid.x> 0 && kid.y> 0);
    
    kid.xs = 0;
    kid.ys = 0;
    
    mutate(kid.x, kid.y);
}

void create_son(Individual &mother, Individual &kid)
{
    assert(mother.x > 0 && mother.y > 0);
    
    kid.x =  mother.x;
    kid.y =  mother.y;
    
    assert(kid.x> 0 && kid.y> 0);
    
    kid.xs = 0;
    kid.ys = 0;
    
    mutate(kid.x, kid.y);
}




//Calculate and register trait values
void displaystats(Population FEMAjuvfree, Population FEMAjuvguarded, Population FEMAsexfree, Population FEMAsexmate, Population MALEguardfree, Population MALEguardpaired, Population MALEsearchfree, Population MALEsexmate)
{
    double xcount=0;
    double ycount=0;
    double stdevsumx=0;
    double stdevsumy=0;
    
    //stats
    for (int i =0; i < NuFjf; i++)
    {
        xcount += FEMAjuvfree[i].x;
        ycount += FEMAjuvfree[i].y;
    }
    for (int i =0; i < NuFjp; i++)
    {
        xcount += FEMAjuvpaired[i].x;
        ycount += FEMAjuvpaired[i].y;
    }
    for (int i =0; i < NuFsf; i++)
    {
        xcount += FEMAsexfree[i].x;
        ycount += FEMAsexfree[i].y;
    }
    for (int i =0; i < NuFsm; i++)
    {
        xcount += FEMAsexmate[i].x;
        ycount += FEMAsexmate[i].y;
    }
    
    for (int i =0; i < NuFmf; i++)
    {
        xcount += FEMAmaturefree[i].x;
        ycount += FEMAmaturefree[i].y;
    }
    for (int i =0; i < NuFmm; i++)
    {
        xcount += FEMAmaturemate[i].x;
        ycount += FEMAmaturemate[i].y;
    }
    
    for (int i =0; i < NuMgf; i++)
    {
        xcount += MALEguardfree[i].x;
        ycount += MALEguardfree[i].y;
    }
    for (int i =0; i < NuMgp; i++)
    {
        xcount += MALEguardpaired[i].x;
        ycount += MALEguardpaired[i].y;
    }
    for (int i =0; i < NuMsf; i++)
    {
        xcount += MALEsearchfree[i].x;
        ycount += MALEsearchfree[i].y;
    }
    for (int i =0; i < NuMsm; i++)
    {
        xcount += MALEsexmate[i].x;
        ycount += MALEsexmate[i].y;
    }
    meanx = xcount / (Npop);
    meany = ycount / (Npop);
    
    for (int i =0; i < NuFjf; i++)
    {
        stdevsumx += pow(FEMAjuvfree[i].x - meanx,2);
        stdevsumy += pow(FEMAjuvfree[i].y - meany,2);
    }
    for (int i =0; i < NuFjp; i++)
    {
        stdevsumx += pow(FEMAjuvpaired[i].x - meanx,2);
        stdevsumy += pow(FEMAjuvpaired[i].y - meany,2);
    }
    for (int i =0; i < NuFsf; i++)
    {
        stdevsumx += pow(FEMAsexfree[i].x - meanx,2);
        stdevsumy += pow(FEMAsexfree[i].y - meany,2);
    }
    for (int i =0; i < NuFsm; i++)
    {
        stdevsumx += pow(FEMAsexmate[i].x - meanx,2);
        stdevsumy += pow(FEMAsexmate[i].y - meany,2);
    }
    
    
    for (int i =0; i < NuFmf; i++)
    {
        stdevsumx += pow(FEMAmaturefree[i].x - meanx,2);
        stdevsumy += pow(FEMAmaturefree[i].y - meany,2);
    }
    for (int i =0; i < NuFmm; i++)
    {
        stdevsumx += pow(FEMAmaturemate[i].x - meanx,2);
        stdevsumy += pow(FEMAmaturemate[i].y - meany,2);
    }
    
    for (int i =0; i < NuMgf; i++)
    {
        stdevsumx += pow(MALEguardfree[i].x - meanx,2);
        stdevsumy += pow(MALEguardfree[i].y - meany,2);
    }
    for (int i =0; i < NuMgp; i++)
    {
        stdevsumx += pow(MALEguardpaired[i].x - meanx,2);
        stdevsumy += pow(MALEguardpaired[i].y - meany,2);
    }
    for (int i =0; i < NuMsf; i++)
    {
        stdevsumx += pow(MALEsearchfree[i].x - meanx,2);
        stdevsumy += pow(MALEsearchfree[i].y - meany,2);
    }
    for (int i =0; i < NuMsm; i++)
    {
        stdevsumx += pow(MALEsexmate[i].x - meanx,2);
        stdevsumy += pow(MALEsexmate[i].y - meany,2);
    }
    
    stdevx=stdevsumx/Npop;
    stdevy=stdevsumy/Npop;
}


//  Main function
int main(int argc, char ** argv)
{
    double meanx_mean,meanx_sum,meany_mean,meany_sum,stdevx_mean,stdevx_sum,stdevy_mean,stdevy_sum;
    int reps;
    double time;
    

    //implement input arguments and intial class values.
    initArguments(argc, argv);
    rng_r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng_r, seed);
    
    FILE *ptt1;
    FILE *ptt2;
    char arq1[memsize];
    char arq2[memsize];
    
    cout <<"simplemites3... start!"<<endl;

    
    meanx_mean=meanx_sum=meany_mean=meany_sum=stdevx_mean=stdevx_sum=stdevy_mean=stdevy_sum=0;
    //starting reps
    for (reps = 0; reps < repmax; reps++)
    {
       
        Init();
        cout <<"rep is "<<reps<<endl;
        time=0;
        events=0;
        int its=0;
        int index=0;

        switcher=1;
        //starting events in each rep
        while (time < tmax )
        {
            do_stats=switcher;

            //define when to do stats
            if(index <= time)
            {   
                switcher = 1;
                cout <<"printing at time "<<time<<endl;
                cout <<reps<<" "<<time<<" "<<"events: "<<events<<"means: " <<meanx<<" "<<stdevx<<" " <<meany<<" "<<stdevy<<endl;
                its++;
                index=skip*its;
                events=0;
            }

                       
            //do_stats = time % skip == 0;
            
            if(NuFjf+NuFjp+NuFsf+NuFsm+NuFmf+NuFmm+NuMgf+NuMgp+NuMsf+NuMsm==0)
            {
                cout <<"Population Extinction!!! at t= "<<time<<endl;
                break;
            }
            
            
            
            assert((NuFjf+NuFjp+NuFsf+NuFsm+NuFmf+NuFmm+NuMgf+NuMgp+NuMsf+NuMsm) <= Npop);
            assert(NuFjf >= 0 && NuFjf < Npop);
            assert(NuFjp >= 0 && NuFjp < Npop);
            assert(NuFsf >= 0 && NuFsf < Npop);
            assert(NuFsm >= 0 && NuFsm < Npop);
            assert(NuFmf >= 0 && NuFmf < Npop);
            assert(NuFmm >= 0 && NuFmm < Npop);
            assert(NuMgf >= 0 && NuMgf < Npop);
            assert(NuMgp >= 0 && NuMgp < Npop);
            assert(NuMsf >= 0 && NuMsf < Npop);
            assert(NuMsm >= 0 && NuMsm < Npop);
            assert(NuMgp == NuFjp);
            
            //assert all classes have "normal" values
            for(int i = 0 ; i < NuFjf; i++) assert(FEMAjuvfree[i].x > 0 && FEMAjuvfree[i].y >0 && FEMAjuvfree[i].xs >= 0 && FEMAjuvfree[i].ys >= 0);
            for(int i = 0 ; i < NuFjp; i++) assert(FEMAjuvpaired[i].x > 0 && FEMAjuvpaired[i].y >0 && FEMAjuvpaired[i].xs >= 0 && FEMAjuvpaired[i].ys >= 0);
            for(int i = 0 ; i < NuFsf; i++) assert(FEMAsexfree[i].x > 0 && FEMAsexfree[i].y >0 && FEMAsexfree[i].xs >= 0 && FEMAsexfree[i].ys >= 0);
            for(int i = 0 ; i < NuFsm; i++) assert(FEMAsexmate[i].x > 0 && FEMAsexmate[i].y >0 && FEMAsexmate[i].xs >= 0 && FEMAsexmate[i].ys >= 0);
            for(int i = 0 ; i < NuFmf; i++) assert(FEMAmaturefree[i].x > 0 && FEMAmaturefree[i].y >0 && FEMAmaturefree[i].xs >= 0 && FEMAmaturefree[i].ys >= 0);
            for(int i = 0 ; i < NuFmm; i++) assert(FEMAmaturemate[i].x > 0 && FEMAmaturemate[i].y >0 && FEMAmaturemate[i].xs >= 0 && FEMAmaturemate[i].ys >= 0);
            for(int i = 0 ; i < NuMgf; i++) assert(MALEguardfree[i].x > 0 && MALEguardfree[i].y > 0);
            for(int i = 0 ; i < NuMgp; i++) assert(MALEguardpaired[i].x > 0 && MALEguardpaired[i].y > 0 );
            for(int i = 0 ; i < NuMsf; i++) assert(MALEsearchfree[i].x > 0 && MALEsearchfree[i].y > 0 );
            for(int i = 0 ; i < NuMsm; i++) assert(MALEsexmate[i].x > 0 && MALEsexmate[i].y > 0 );
            
            //determine the number of possible events
            nprobs = 22;
            //define vector of possible events
            double *probs = new double [nprobs];
            //define sum of all probabilities for cummulative distribution
            double sumprobs=0;
            
            //DETERMINE PROBS FOR ALL EVENTS:
            //1. maturity of uFjf
            sumprobs = probs[0] = NuFjf*(1/Tf);//probability of uFjf->uFsf
            assert(sumprobs>=0);
            if(NuFjf==0) assert(probs[0]==0);
            
            //2. maturity of uFjp
            sumprobs = probs[1] = sumprobs + NuFjp*(1/Tf);//probability of uFjp->uFmf
            assert(sumprobs>=0);
            if(NuFjp==0) assert(probs[1]==probs[0]);
            
            //3.abandon guarding
            sumprobs = probs[2] = sumprobs + NuMgp * (1/Tg);
            assert(sumprobs>=0);
            if(NuMgp==0) assert(probs[2]==probs[1]);

            //4. sexual rematuration of uFsm
            sumprobs = probs[3] = sumprobs + NuFsm*(1 / Tfo);//probability of uFmf->uFsf
            assert(sumprobs>=0);
            if(NuFsm==0) assert(probs[3]==probs[2]);
            
            //5. sexual rematuration of uFmm
            sumprobs = probs[4] = sumprobs + NuFmm*(1 / Tfo);//probability of uFmm->uFmf
            assert(sumprobs>=0);
            if(NuFmm==0) assert(probs[4]==probs[3]);
            
            //6. sexual rematuration of uMsm
            sumprobs = probs[5] = sumprobs + NuMsm*(1 / Tmo);//probability of uMmf->uMsf
            assert(sumprobs>=0);
            if(NuMsm==0) assert(probs[5]==probs[4]);
            
            //7. death of uFjf
            sumprobs = probs[6] = sumprobs + NuFjf*muFjf;
            assert(sumprobs>=0);
            if(NuFjf==0) assert(probs[6]==probs[5]);
            
            
            //8. death of uFjp
            sumprobs = probs[7] = sumprobs + NuFjp*muFjp;
            assert(sumprobs>=0);
            if(NuFjp==0) assert(probs[7]==probs[6]);
            
            
            //9. death of uFsf
            sumprobs = probs[8] = sumprobs + NuFsf*muFsf;
            assert(sumprobs>=0);
            if(NuFsf==0) assert(probs[8]==probs[7]);
            
            
            //10. death of uFsm
            sumprobs = probs[9] = sumprobs + NuFsm*muFsm;
            assert(sumprobs>=0);
            if(NuFsm==0) assert(probs[9]==probs[8]);
            
            
            //11. death of uFmf
            sumprobs = probs[10] = sumprobs + NuFmf*muFmf;
            assert(sumprobs>=0);
            if(NuFmf==0) assert(probs[10]==probs[9]);
            
            //12. death of uFmm
            sumprobs = probs[11] = sumprobs + NuFmm*muFmm;
            assert(sumprobs>=0);
            if(NuFmm==0) assert(probs[11]==probs[10]);
            
            
            //13. death of uMgf
            sumprobs = probs[12] = sumprobs + NuMgf*muMgf;
            assert(sumprobs>=0);
            if(NuMgf==0) assert(probs[12]==probs[11]);
            
            
            //14. death of uMgp
            sumprobs = probs[13] = sumprobs + NuMgp*muMgp*cg;
            assert(sumprobs>=0);
            if(NuMgp==0) assert(probs[13]==probs[12]);
            
            
            //15. death of uMsf
            sumprobs = probs[14] = sumprobs + NuMsf*muMsf*cs;
            assert(sumprobs>=0);
            if(NuMsf==0) assert(probs[14]==probs[13]);
            
            //            cout << "sumprobs bf " <<sumprobs<<"\t";
            
            //16. death of uMsm
            sumprobs = probs[15] = sumprobs + NuMsm*muMsm;
            assert(sumprobs>=0);
            if(NuMsm==0) assert(probs[15]==probs[14]);
            //            cout << "sumprobs af " <<sumprobs<<endl;
            
            //17. male guarder finds juvenile female uMgf->uFjf
            if(NuMgf>0){
                sumprobs = probs[16] = sumprobs + NuMgf*sqrt(NuFjf/NuMgf);
            }
            else{
                    sumprobs = probs[16] = sumprobs;
            }    
            assert(sumprobs>=0);
            if(NuMgf==0) assert(probs[16]==probs[15]);

            
            //18. juvenile female finds guarder male
            if(NuFjf>0){
                sumprobs = probs[17] = sumprobs + NuFjf*sqrt(NuMgf/NuFjf);
            }
            else{
                sumprobs = probs[17] = sumprobs;
            }
            assert(sumprobs>=0);
            if(NuFjf==0) assert(probs[17]==probs[16]);

            //19. searcher male finds virgin female
            if(NuMsf>0){
                sumprobs = probs[18] = sumprobs + NuMsf*sqrt(NuFsf/NuMsf);
            }
            else{
                sumprobs = probs[18] = sumprobs;
            }
            assert(sumprobs>=0);
            if(NuMsf==0) assert(probs[18]==probs[17]);

            
            //20. virgin female finds searcher male
            if(NuFsf>0){
                sumprobs = probs[19] = sumprobs + NuFsf*sqrt(NuMsf/NuFsf);
            }
            else{
                sumprobs = probs[19] = sumprobs;
            }
            assert(sumprobs>=0);
            if(NuFsf==0) assert(probs[19]==probs[18]);


            //21. searcher male finds virgin female
            if(NuMsf>0){
                sumprobs = probs[20] = sumprobs + NuMsf*sqrt(NuFmf/NuMsf);
            }
            else{
                sumprobs = probs[20] = sumprobs;
            }
            assert(sumprobs>=0);
            if(NuMsf==0) assert(probs[20]==probs[19]);

            
            //22. virgin female finds searcher male
            if(NuFmf>0){
                sumprobs = probs[21] = sumprobs + NuFmf*sqrt(NuMsf/NuFmf);
            }
            else{
                sumprobs = probs[21] = sumprobs;
            }
            assert(sumprobs>=0);
            if(NuFmf==0) assert(probs[21]==probs[20]);

            sum_events=22;
            
            //DETERMINE WINNING EVENT:
            //draw a probability out of the cumulative probability distribution
            
            //            double randvalue = gsl_rng_uniform(rng_r);
            //            prob = randvalue * sumprobs;
            double prob = gsl_ran_flat(rng_r, 0, sumprobs);
            
            // choose the value
            int choice=0;
            for (choice = 0; choice < nprobs; choice++)
            {
                assert(probs[choice]<=sumprobs);
                // value found, break
                if (prob <= probs[choice])
                {
                    break;
                }
            }
            
            
            //delete array of probabilities (takes too much memory)
            delete[] probs;
            
            
            //Now, we need to find to what [choice] corresponds to.
            
            //1. maturity of uFjf. When they mature they may produce offspring
            if (choice < 1)
            {
                int indiv = int(gsl_rng_uniform(rng_r) * NuFjf);
                
                FEMAsexfree[NuFsf].x = FEMAjuvfree[indiv].x;
                FEMAsexfree[NuFsf].y = FEMAjuvfree[indiv].y;
                FEMAsexfree[NuFsf].xs = FEMAjuvfree[indiv].xs;
                FEMAsexfree[NuFsf].ys = FEMAjuvfree[indiv].ys;
                assert(FEMAsexfree[NuFsf].x > 0 && FEMAsexfree[NuFsf].y > 0);
                assert(FEMAsexfree[NuFsf].x > 0 || FEMAsexfree[NuFsf].y > 0);
                NuFsf++;
                
                FEMAjuvfree[indiv].x = FEMAjuvfree[NuFjf - 1].x;
                FEMAjuvfree[indiv].y = FEMAjuvfree[NuFjf - 1].y;
                FEMAjuvfree[indiv].xs = FEMAjuvfree[NuFjf - 1].xs;
                FEMAjuvfree[indiv].ys = FEMAjuvfree[NuFjf - 1].ys;
                assert(FEMAjuvfree[indiv].x > 0 && FEMAjuvfree[indiv].y > 0);
                NuFjf--;
                
                //Assexual reproduction, where female produces as many offspring as necessary to get back to Npop
                
                if (gsl_rng_uniform(rng_r) < assex)
                {
                    for (int newborn_i = 0; newborn_i < Npop - (NuFjf + NuFjp + NuFsf + NuFsm + NuFmf + NuFmm + NuMgf + NuMgp + NuMsf + NuMsm); ++newborn_i)
                    {
                        Individual kid; //selection of strategy
                        create_son(FEMAsexfree[NuFsf-1], kid);
                        if (gsl_rng_uniform(rng_r) < kid.x)
                        {
                            MALEguardfree[NuMgf++]=kid;
                        }
                        else
                        {
                            MALEsearchfree[NuMsf++]=kid;
                        }
                    }
                }
            }
            
            //2. maturity of uFjp >> mating of NuFjp and NuMgp
            else if (choice < 2)
            {
                
                int indiv = int(gsl_rng_uniform(rng_r) * NuFjp);
                
                //choose male guarding to mate with female.
                int pair = int(gsl_rng_uniform(rng_r) * NuMgp);
                
                //update female classes
                FEMAsexmate[NuFsm].x = FEMAjuvpaired[indiv].x;
                FEMAsexmate[NuFsm].y = FEMAjuvpaired[indiv].y;
                FEMAsexmate[NuFsm].xs = MALEguardpaired[pair].x;
                FEMAsexmate[NuFsm].ys = MALEguardpaired[pair].y;
                assert(FEMAsexmate[NuFsm].x>0 && FEMAsexmate[NuFsm].y > 0 && FEMAsexmate[NuFsm].xs > 0 && FEMAsexmate[NuFsm].ys > 0);
                NuFsm++;
                
                FEMAjuvpaired[indiv].x = FEMAjuvpaired[NuFjp - 1].x;
                FEMAjuvpaired[indiv].y = FEMAjuvpaired[NuFjp - 1].y;
                FEMAjuvpaired[indiv].xs = FEMAjuvpaired[NuFjp - 1].xs;
                FEMAjuvpaired[indiv].ys = FEMAjuvpaired[NuFjp - 1].ys;
                assert(FEMAjuvpaired[indiv].x>0 && FEMAjuvpaired[indiv].y>0);
                NuFjp--;
                
                //update male classes
                MALEsexmate[NuMsm].x = MALEguardpaired[pair].x;
                MALEsexmate[NuMsm].y = MALEguardpaired[pair].y;
                MALEsexmate[NuMsm].xs = MALEguardpaired[pair].xs;
                MALEsexmate[NuMsm].ys = MALEguardpaired[pair].ys;
                assert(MALEsexmate[NuMsm].x>0 && MALEsexmate[NuMsm].y>0);
                NuMsm++;
                
                MALEguardpaired[pair].x = MALEguardpaired[NuMgp - 1].x;
                MALEguardpaired[pair].y = MALEguardpaired[NuMgp - 1].y;
                MALEguardpaired[pair].xs = MALEguardpaired[NuMgp - 1].xs;
                MALEguardpaired[pair].ys = MALEguardpaired[NuMgp - 1].ys;
                assert(MALEguardpaired[pair].x>0 && MALEguardpaired[pair].y>0);
                NuMgp--;
                
                // pair produces as many offspring as necessary to get back to Npop
                for (int newborn_i = 0;
                     newborn_i < Npop - (NuFjf + NuFjp + NuFsf + NuFsm + NuFmf + NuFmm + NuMgf + NuMgp + NuMsf + NuMsm);
                     ++newborn_i)
                {
                    
                    Individual kid;
                    create_kid(FEMAsexmate[NuFsm-1], kid);
                    
                    //assign kid to either male or female "IN" group
                    if (gsl_rng_uniform(rng_r) > sexratio1)
                    {
                        FEMAjuvfree[NuFjf++] = kid;
                    }
                    else
                    {
                        //assign male kid to either searching or guarding strategies
                        if (gsl_rng_uniform(rng_r) < kid.x)
                        {
                            MALEguardfree[NuMgf++]=kid;
                        }
                        else
                        {
                            MALEsearchfree[NuMsf++]=kid;
                        }
                    }
                    //                    cout <<"912 newborn from mature of guarding female"<<endl;
                    //                    cout <<"912 conditions after reproduction: "<<Npop<<" "<<NuMgf + NuMgp + NuMsf + NuMsm + NuFjf + NuFjp + NuFsf + NuFsm<<endl;
                    
                }
                
            }
            
            //3.abandon guarding
            else if (choice < 3)
            {
                int indiv = int(gsl_rng_uniform(rng_r) * NuMgp);
                
                //choose female pair
                int pair = int(gsl_rng_uniform(rng_r)*NuFjp);
                assert(pair < NuFjp);
                
                FEMAjuvfree[NuFjf].x = FEMAjuvpaired[pair].x;
                FEMAjuvfree[NuFjf].y = FEMAjuvpaired[pair].y;
                FEMAjuvfree[NuFjf].xs = FEMAjuvpaired[pair].xs;
                FEMAjuvfree[NuFjf].ys = FEMAjuvpaired[pair].ys;
                assert(FEMAjuvfree[NuFjf].x > 0 && FEMAjuvfree[NuFjf].y > 0);
                NuFjf++;
                
                FEMAjuvpaired[pair].x = FEMAjuvpaired[NuFjp - 1].x;
                FEMAjuvpaired[pair].y = FEMAjuvpaired[NuFjp - 1].y;
                FEMAjuvpaired[pair].xs = FEMAjuvpaired[NuFjp - 1].xs;
                FEMAjuvpaired[pair].ys = FEMAjuvpaired[NuFjp - 1].ys;
                assert(FEMAjuvpaired[pair].x > 0 && FEMAjuvpaired[pair].y > 0);
                NuFjp--;
                
                //update male free class
                MALEguardfree[NuMgf].x = MALEguardpaired[indiv].x;
                MALEguardfree[NuMgf].y = MALEguardpaired[indiv].y;
                MALEguardfree[NuMgf].xs = MALEguardpaired[indiv].xs;
                MALEguardfree[NuMgf].ys = MALEguardpaired[indiv].ys;
                assert(MALEguardfree[NuMgf].x > 0 && MALEguardfree[NuMgf].y > 0);
                NuMgf++;
                
                //update male guarding class
                MALEguardpaired[indiv].x = MALEguardpaired[NuMgp - 1].x;
                MALEguardpaired[indiv].y = MALEguardpaired[NuMgp - 1].y;
                MALEguardpaired[indiv].xs = MALEguardpaired[NuMgp - 1].xs;
                MALEguardpaired[indiv].ys = MALEguardpaired[NuMgp - 1].ys;
                assert(MALEguardpaired[indiv].x > 0 && MALEguardpaired[indiv].y > 0);
                NuMgp--;

            }
            
            
            
            
            
            //4. sexual re-maturation of uFsm
            else if (choice < 4)
            {
                int indiv = int(gsl_rng_uniform(rng_r) * NuFsm);
                
                //update female classes
                FEMAmaturefree[NuFmf].x = FEMAsexmate[indiv].x;
                FEMAmaturefree[NuFmf].y = FEMAsexmate[indiv].y;
                FEMAmaturefree[NuFmf].xs = FEMAsexmate[indiv].xs;
                FEMAmaturefree[NuFmf].ys = FEMAsexmate[indiv].ys;
                assert(FEMAmaturefree[NuFmf].x>0 && FEMAmaturefree[NuFmf].y>0 && FEMAmaturefree[NuFmf].xs>0 && FEMAmaturefree[NuFmf].ys>0);
                NuFmf++;
                
                FEMAsexmate[indiv].x = FEMAsexmate[NuFsm - 1].x;
                FEMAsexmate[indiv].y = FEMAsexmate[NuFsm - 1].y;
                FEMAsexmate[indiv].xs = FEMAsexmate[NuFsm - 1].xs;
                FEMAsexmate[indiv].ys = FEMAsexmate[NuFsm - 1].ys;
                assert(FEMAsexmate[indiv].x>0 && FEMAsexmate[indiv].y>0);
                NuFsm--;
            }
            
            //5. sexual re-maturation of uFmm
            else if (choice < 5)
            {
                int indiv = int(gsl_rng_uniform(rng_r) * NuFmm);
                
                //update female classes
                FEMAmaturefree[NuFmf].x = FEMAmaturemate[indiv].x;
                FEMAmaturefree[NuFmf].y = FEMAmaturemate[indiv].y;
                FEMAmaturefree[NuFmf].xs = FEMAmaturemate[indiv].xs;
                FEMAmaturefree[NuFmf].ys = FEMAmaturemate[indiv].ys;
                assert(FEMAmaturefree[NuFmf].x>0 && FEMAmaturefree[NuFmf].y>0);
                NuFmf++;
                
                FEMAmaturemate[indiv].x = FEMAmaturemate[NuFmm - 1].x;
                FEMAmaturemate[indiv].y = FEMAmaturemate[NuFmm - 1].y;
                FEMAmaturemate[indiv].xs = FEMAmaturemate[NuFmm - 1].xs;
                FEMAmaturemate[indiv].ys = FEMAmaturemate[NuFmm - 1].ys;
                assert(FEMAmaturemate[indiv].x>0 && FEMAmaturemate[indiv].y>0);
                NuFmm--;
            }
            
            //6. sexual rematuration of uMsm
            // This is the main reason why we are limiting x to vary between 0 and 1. Male function is determined by drawing a random number within these limits and comparing it to the value of x.
            else if (choice < 6)
            {
                int indiv = int(gsl_rng_uniform(rng_r) * NuMsm);
                
                //Now determine if this male will be searcher or guarder
                if(gsl_rng_uniform(rng_r) < MALEsexmate[indiv].x)                  {
                    //male is guarder
                    
                    //update male classes
                    MALEguardfree[NuMgf].x = MALEsexmate[indiv].x;
                    MALEguardfree[NuMgf].y = MALEsexmate[indiv].y;
                    MALEguardfree[NuMgf].xs = MALEsexmate[indiv].xs;
                    MALEguardfree[NuMgf].ys = MALEsexmate[indiv].ys;
                    assert(MALEguardfree[NuMgf].x>0 && MALEguardfree[NuMgf].y>0);
                    NuMgf++;
                    
                    MALEsexmate[indiv].x = MALEsexmate[NuMsm-1].x;
                    MALEsexmate[indiv].y = MALEsexmate[NuMsm-1].y;
                    MALEsexmate[indiv].xs = MALEsexmate[NuMsm-1].xs;
                    MALEsexmate[indiv].ys = MALEsexmate[NuMsm-1].ys;
                    assert(MALEsexmate[indiv].x>0 && MALEsexmate[indiv].y>0);
                    NuMsm--;
                }
                else
                {
                    //male is searcher
                    
                    //update male classes
                    MALEsearchfree[NuMsf].x = MALEsexmate[indiv].x;
                    MALEsearchfree[NuMsf].y = MALEsexmate[indiv].y;
                    MALEsearchfree[NuMsf].xs = MALEsexmate[indiv].xs;
                    MALEsearchfree[NuMsf].ys = MALEsexmate[indiv].ys;
                    assert(MALEsearchfree[NuMsf].x>0 && MALEsearchfree[NuMsf].y>0);
                    NuMsf++;
                    
                    MALEsexmate[indiv].x = MALEsexmate[NuMsm - 1].x;
                    MALEsexmate[indiv].y = MALEsexmate[NuMsm - 1].y;
                    MALEsexmate[indiv].xs = MALEsexmate[NuMsm - 1].xs;
                    MALEsexmate[indiv].ys = MALEsexmate[NuMsm - 1].ys;
                    assert(MALEsexmate[indiv].x>0 && MALEsexmate[indiv].y>0);
                    NuMsm--;
                }
            }
            
            //7. death of uFjf
            else if (choice < 7)
            {
                int indiv = int(gsl_rng_uniform(rng_r)*NuFjf);
                int newvalue;
                newvalue = death_FEMAjuvfree(FEMAjuvfree, indiv, NuFjf);
                NuFjf=newvalue;
            }
            
            //8. death of uFjp
            else if (choice < 8)
            {
                int indiv = int(gsl_rng_uniform(rng_r)*NuFjp);
                
                //choose male pair
                int pair = int(gsl_rng_uniform(rng_r)*NuMgp);
                
                MALEguardfree[NuMgf].x = MALEguardpaired[pair].x;
                MALEguardfree[NuMgf].y = MALEguardpaired[pair].y;
                MALEguardfree[NuMgf].xs = MALEguardpaired[pair].xs;
                MALEguardfree[NuMgf].ys = MALEguardpaired[pair].ys;
                assert(MALEguardfree[NuMgf].x > 0 && MALEguardfree[NuMgf].y > 0);
                assert(MALEguardfree[NuMgf].x == MALEguardpaired[pair].x && MALEguardfree[NuMgf].y ==MALEguardpaired[pair].y);
                NuMgf++;
                
                MALEguardpaired[pair].x = MALEguardpaired[NuMgp - 1].x;
                MALEguardpaired[pair].y = MALEguardpaired[NuMgp - 1].y;
                MALEguardpaired[pair].xs = MALEguardpaired[NuMgp - 1].xs;
                MALEguardpaired[pair].ys = MALEguardpaired[NuMgp - 1].ys;
                assert(MALEguardpaired[pair].x > 0 && MALEguardpaired[pair].y > 0);
                NuMgp--;
                
                int newvalue;
                newvalue = death_FEMAjuvpaired(FEMAjuvpaired, indiv, NuFjp);
                NuFjp=newvalue;
            }
            
            //9. death of uFsf
            else if (choice < 9)
            {
                int indiv = int(gsl_rng_uniform(rng_r)*NuFsf);
                
                int newvalue;
                newvalue = death_FEMAsexfree(FEMAsexfree, indiv, NuFsf);
                NuFsf=newvalue;
            }
            
            //10. death of uFsm
            else if (choice < 10)
            {
                int indiv = int(gsl_rng_uniform(rng_r)*NuFsm);
                
                int newvalue;
                newvalue = death_FEMAsexmate(FEMAsexmate, indiv, NuFsm);
                NuFsm = newvalue;
            }
            
            //11. death of uFmf
            else if (choice < 11)
            {
                int indiv = int(gsl_rng_uniform(rng_r)*NuFmf);
                
                int newvalue;
                newvalue = death_FEMAmaturefree(FEMAmaturefree, indiv, NuFmf);
                NuFmf=newvalue;
            }
            
            //12. death of uFmm
            else if (choice < 12)
            {
                int indiv = int(gsl_rng_uniform(rng_r)*NuFmm);
                
                int newvalue;
                newvalue = death_FEMAmaturemate(FEMAmaturemate, indiv, NuFmm);
                NuFmm=newvalue;
            }
            
            //13. death of uMgf
            else if (choice < 13)
            {
                int indiv = int(gsl_rng_uniform(rng_r)*NuMgf);
                
                int newvalue;
                newvalue = death_MALEguardfree(MALEguardfree, indiv, NuMgf);
                NuMgf=newvalue;
                
            }
            
            //14. death of uMgp
            else if (choice < 14)
            {
                int indiv = int(gsl_rng_uniform(rng_r)*NuMgp);
                
                //choose female pair
                int pair = int(gsl_rng_uniform(rng_r)*NuFjp);
                assert(pair < NuFjp);
                
                FEMAjuvfree[NuFjf].x = FEMAjuvpaired[pair].x;
                FEMAjuvfree[NuFjf].y = FEMAjuvpaired[pair].y;
                FEMAjuvfree[NuFjf].xs = FEMAjuvpaired[pair].xs;
                FEMAjuvfree[NuFjf].ys = FEMAjuvpaired[pair].ys;
                assert(FEMAjuvfree[NuFjf].x > 0 && FEMAjuvfree[NuFjf].y > 0);
                NuFjf++;
                
                FEMAjuvpaired[pair].x = FEMAjuvpaired[NuFjp - 1].x;
                FEMAjuvpaired[pair].y = FEMAjuvpaired[NuFjp - 1].y;
                FEMAjuvpaired[pair].xs = FEMAjuvpaired[NuFjp - 1].xs;
                FEMAjuvpaired[pair].ys = FEMAjuvpaired[NuFjp - 1].ys;
                assert(FEMAjuvpaired[pair].x > 0 && FEMAjuvpaired[pair].y > 0);
                NuFjp--;
                
                int newvalue;
                newvalue = death_MALEguardpaired(MALEguardpaired, indiv, NuMgp);
                NuMgp=newvalue;
            }
            
            
            //15. death of uMsf
            else if (choice < 15)
            {
                int indiv = int(gsl_rng_uniform(rng_r)*NuMsf);
                
                int newvalue;
                newvalue = death_MALEsearchfree(MALEsearchfree, indiv, NuMsf);
                NuMsf=newvalue;
            }
            
            //16. death of uMsm
            else if (choice < 16)
            {
                int indiv = int(gsl_rng_uniform(rng_r)*NuMsm);
                //                cout<<"choice "<<choice<<endl;
                int newvalue;
                newvalue = death_MALEsexmate(MALEsexmate, indiv, NuMsm);
                NuMsm=newvalue;
            }
            
            //17&18. Pairing of uFjf/uMgf
            else if (choice < 18)
            {
                indivF = int(gsl_rng_uniform(rng_r)*NuFjf);
                indivM = int(gsl_rng_uniform(rng_r)*NuMgf);
                assert(indivF < NuFjf);
                assert(indivM < NuMgf);
                assert(FEMAjuvfree[indivF].x>0 && FEMAjuvfree[indivF].y>0);
                assert(MALEguardfree[indivM].x>0 && MALEguardfree[indivM].y>0);
                
                //female receives male sperm
                FEMAjuvfree[indivF].xs = MALEguardfree[indivM].x;
                FEMAjuvfree[indivF].ys = MALEguardfree[indivM].y;
                
                assert(FEMAjuvfree[indivF].x > 0 && FEMAjuvfree[indivF].y > 0);
                assert(MALEguardfree[indivM].x > 0 && MALEguardfree[indivM].y > 0);
                assert(FEMAjuvfree[indivF].x > 0 && FEMAjuvfree[indivF].y > 0 && FEMAjuvfree[indivF].xs > 0 && FEMAjuvfree[indivF].ys > 0);
                
                // update guarding classes
                FEMAjuvpaired[NuFjp].x = FEMAjuvfree[indivF].x;
                FEMAjuvpaired[NuFjp].y = FEMAjuvfree[indivF].y;
                FEMAjuvpaired[NuFjp].xs = FEMAjuvfree[indivF].xs;
                FEMAjuvpaired[NuFjp].ys = FEMAjuvfree[indivF].ys;
                assert(FEMAjuvpaired[NuFjp].x == FEMAjuvfree[indivF].x);
                assert(FEMAjuvpaired[NuFjp].y == FEMAjuvfree[indivF].y);
                assert(FEMAjuvpaired[NuFjp].x > 0 && FEMAjuvpaired[NuFjp].y > 0 && FEMAjuvpaired[NuFjp].xs > 0 && FEMAjuvpaired[NuFjp].ys > 0);
                NuFjp++;
                
                MALEguardpaired[NuMgp].x = MALEguardfree[indivM].x;
                MALEguardpaired[NuMgp].y = MALEguardfree[indivM].y;
                MALEguardpaired[NuMgp].xs = MALEguardfree[indivM].xs;
                MALEguardpaired[NuMgp].ys = MALEguardfree[indivM].ys;
                assert(MALEguardpaired[NuMgp].x == MALEguardfree[indivM].x);
                assert(MALEguardpaired[NuMgp].y == MALEguardfree[indivM].y);
                assert(MALEguardpaired[NuMgp].x > 0 && MALEguardpaired[NuMgp].y > 0);
                NuMgp++;
                
                //update free classes
                FEMAjuvfree[indivF].x = FEMAjuvfree[NuFjf - 1].x;
                FEMAjuvfree[indivF].y = FEMAjuvfree[NuFjf - 1].y;
                FEMAjuvfree[indivF].xs = FEMAjuvfree[NuFjf - 1].xs;
                FEMAjuvfree[indivF].ys = FEMAjuvfree[NuFjf - 1].ys;
                assert(FEMAjuvfree[indivF].x > 0 && FEMAjuvfree[indivF].y > 0);
                NuFjf--;
                
                MALEguardfree[indivM].x = MALEguardfree[NuMgf - 1].x;
                MALEguardfree[indivM].y = MALEguardfree[NuMgf - 1].y;
                MALEguardfree[indivM].xs = MALEguardfree[NuMgf - 1].xs;
                MALEguardfree[indivM].ys = MALEguardfree[NuMgf - 1].ys;
                assert(MALEguardfree[indivM].x > 0 && MALEguardfree[indivM].y > 0);
                NuMgf--;
            }
            
            //19&20. Pairing and mating of uFsf/uMsf
            else if (choice < 20)
            {
                indivF = int(gsl_rng_uniform(rng_r)*NuFsf);
                indivM = int(gsl_rng_uniform(rng_r)*NuMsf);

                assert(indivF < NuFsf);
                assert(indivM < NuMsf);
                
                //female receives male sperm
                FEMAsexfree[indivF].xs = MALEsearchfree[indivM].x;
                FEMAsexfree[indivF].ys = MALEsearchfree[indivM].y;
                assert(FEMAsexfree[indivF].x > 0 && FEMAsexfree[indivF].y > 0 && FEMAsexfree[indivF].xs > 0 && FEMAsexfree[indivF].ys > 0);
                
                // update female classes
                FEMAsexmate[NuFsm].x = FEMAsexfree[indivF].x;
                FEMAsexmate[NuFsm].y = FEMAsexfree[indivF].y;
                FEMAsexmate[NuFsm].xs = FEMAsexfree[indivF].xs;
                FEMAsexmate[NuFsm].ys = FEMAsexfree[indivF].ys;
                assert(FEMAsexmate[NuFsm].x > 0 && FEMAsexmate[NuFsm].y > 0 && FEMAsexmate[NuFsm].xs >= 0 && FEMAsexmate[NuFsm].ys>= 0);
                NuFsm++;
                
                FEMAsexfree[indivF].x = FEMAsexfree[NuFsf - 1].x;
                FEMAsexfree[indivF].y = FEMAsexfree[NuFsf - 1].y;
                FEMAsexfree[indivF].xs = FEMAsexfree[NuFsf - 1].xs;
                FEMAsexfree[indivF].ys = FEMAsexfree[NuFsf - 1].ys;
                assert(FEMAsexfree[indivF].x > 0 && FEMAsexfree[indivF].y > 0);
                NuFsf--;
                
                // update male classes
                MALEsexmate[NuMsm].x = MALEsearchfree[indivM].x;
                MALEsexmate[NuMsm].y = MALEsearchfree[indivM].y;
                MALEsexmate[NuMsm].xs = MALEsearchfree[indivM].xs;
                MALEsexmate[NuMsm].ys = MALEsearchfree[indivM].ys;
                assert(MALEsexmate[NuMsm].x > 0 && MALEsexmate[NuMsm].y > 0);
                NuMsm++;
                
                MALEsearchfree[indivM].x = MALEsearchfree[NuMsf - 1].x;
                MALEsearchfree[indivM].y = MALEsearchfree[NuMsf - 1].y;
                MALEsearchfree[indivM].xs = MALEsearchfree[NuMsf - 1].xs;
                MALEsearchfree[indivM].ys = MALEsearchfree[NuMsf - 1].ys;
                assert(MALEsearchfree[indivM].x > 0 && MALEsearchfree[indivM].y > 0);
                NuMsf--;
                
                for (int newborn_i = 0;
                     newborn_i < Npop - (NuFjf + NuFjp + NuFsf + NuFsm + NuFmf + NuFmm + NuMgf + NuMgp + NuMsf + NuMsm);
                     ++newborn_i)
                {
                    Individual kid;
                    create_kid(FEMAsexmate[NuFsm-1], kid);
                    
                    //assign kid to either male or female "IN" group
                    if (gsl_rng_uniform(rng_r) > sexratio1)
                    {
                        FEMAjuvfree[NuFjf++] = kid;
                    }
                    else
                    {
                        //assign male kid to either searching or guarding strategies
                        if (gsl_rng_uniform(rng_r) < kid.x)
                        {
                            MALEguardfree[NuMgf++] = kid;
                        }
                        else
                        {
                            MALEsearchfree[NuMsf++]=kid;
                        }
                    }
                }
                
            }
            
            else //pairing of uFmf with uMsf, but no production of offspring
            {
//                cout <<"something is wrong"<<endl;
//                break;
                indivF = int(gsl_rng_uniform(rng_r)*NuFmf);
                indivM = int(gsl_rng_uniform(rng_r)*NuMsf);
        
                assert(indivF < NuFmf);
                assert(indivM < NuMsf);
                
                if(gsl_rng_uniform(rng_r) > discr)
                {
                    //female receives male sperm
                    FEMAmaturefree[indivF].xs = MALEsearchfree[indivM].x;
                    FEMAmaturefree[indivF].ys = MALEsearchfree[indivM].y;
                    assert(FEMAmaturefree[indivF].x > 0 && FEMAmaturefree[indivF].y > 0 && FEMAmaturefree[indivF].xs > 0 && FEMAmaturefree[indivF].ys > 0);
                    
                    // update female classes
                    FEMAmaturemate[NuFmm].x = FEMAmaturefree[indivF].x;
                    FEMAmaturemate[NuFmm].y = FEMAmaturefree[indivF].y;
                    FEMAmaturemate[NuFmm].xs = FEMAmaturefree[indivF].xs;
                    FEMAmaturemate[NuFmm].ys = FEMAmaturefree[indivF].ys;
                    assert(FEMAmaturemate[NuFmm].x > 0 && FEMAmaturemate[NuFmm].y > 0 && FEMAmaturemate[NuFmm].xs >= 0 && FEMAmaturemate[NuFmm].ys >= 0);
                    NuFmm++;
                    
                    FEMAmaturefree[indivF].x = FEMAmaturefree[NuFmf - 1].x;
                    FEMAmaturefree[indivF].y = FEMAmaturefree[NuFmf - 1].y;
                    FEMAmaturefree[indivF].xs = FEMAmaturefree[NuFmf - 1].xs;
                    FEMAmaturefree[indivF].ys = FEMAmaturefree[NuFmf - 1].ys;
                    assert(FEMAmaturefree[indivF].x > 0 && FEMAmaturefree[indivF].y > 0 && FEMAmaturefree[indivF].xs >= 0 && FEMAmaturefree[indivF].ys >= 0);
                    NuFmf--;
                    
                    // update male classes
                    MALEsexmate[NuMsm].x = MALEsearchfree[indivM].x;
                    MALEsexmate[NuMsm].y = MALEsearchfree[indivM].y;
                    MALEsexmate[NuMsm].xs = MALEsearchfree[indivM].xs;
                    MALEsexmate[NuMsm].ys = MALEsearchfree[indivM].ys;
                    assert(MALEsexmate[NuMsm].x > 0 && MALEsexmate[NuMsm].y > 0);
                    NuMsm++;
                    
                    MALEsearchfree[indivM].x = MALEsearchfree[NuMsf - 1].x;
                    MALEsearchfree[indivM].y = MALEsearchfree[NuMsf - 1].y;
                    MALEsearchfree[indivM].xs = MALEsearchfree[NuMsf - 1].xs;
                    MALEsearchfree[indivM].ys = MALEsearchfree[NuMsf - 1].ys;
                    assert(MALEsearchfree[indivM].x > 0 && MALEsearchfree[indivM].y > 0);
                    NuMsf--;
                }
//                else cout<<"at t="<<time<<", male does not mate with mated female"<<endl;
            }
            
             //DO STATS
            if (do_stats)
            {
                displaystats(FEMAjuvfree, FEMAjuvpaired, FEMAsexfree, FEMAsexmate, MALEguardfree, MALEguardpaired, MALEsearchfree, MALEsexmate);

                sprintf(arq2,"dynamics.csv");
                ptt2 = fopen(arq2,"a");
                fprintf(ptt2,"%d\t %f\t %g\t %g\t %g\t %g\t %d\t %d\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\n",reps, time, meanx, stdevx, meany, stdevy, repmax, tmax, mu, stdev, mortall, cs, cg, Ms, Mg, Tf, Tg,Tfo, Tmo, sexratio1, discr, assex, prop,NuFjf,NuFjp,NuFsf,NuFsm,NuFmf,NuFmm,NuMgf,NuMgp,NuMsf,NuMsm);
                fclose(ptt2);
                cout << "time "<<time<<", meanx: "<<meanx<<", meany: "<<meany<<endl;

                switcher = 0;
            }
            
            if(NuFjf+NuFjp+NuFsf+NuFsm+NuFmf+NuFmm+NuMgf+NuMgp+NuMsf+NuMsm==0)
            {
                cout <<"Population Extinction!!! at t= "<<time<<endl;
                break;
            }
            
            time = time + exp(-1/sum_events)/sum_events;
            events++;
            
        } //end of time
        cout<<time<<endl;
        
        displaystats(FEMAjuvfree, FEMAjuvpaired, FEMAsexfree, FEMAsexmate, MALEguardfree, MALEguardpaired, MALEsearchfree, MALEsexmate);

        meanx_sum += meanx;
        meany_sum += meany;
        stdevx_sum += stdevx;
        stdevy_sum += stdevy;
    } //end of reps

    meanx_mean = double(meanx_sum/reps);
    meany_mean = double(meany_sum/reps);
    stdevx_mean = double(stdevx_sum/reps);
    stdevy_mean = double(stdevy_sum/reps);
    
    cout <<"making of summaryheats.csv" <<endl;

    sprintf(arq1,"summaryheats.csv");
    ptt1 = fopen(arq1,"a");
    fprintf(ptt1,"%d\t %f\t %g\t %g\t %g\t %g\t %d\t %d\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n",reps, time, meanx_mean, stdevx_mean, meany_mean, stdevy_mean, repmax, tmax, mu, stdev, mortall, cs, cg, k, Ms, Mg, Tf, Tg,Tfo, Tmo, sexratio1, discr, assex, prop);
    fclose(ptt1);
   

    
    cout <<"out of here" <<endl;
    
}
