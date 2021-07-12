//
// Created by rasa on 06/07/2021.
// Spatial SIR model
//
#include "Aboria.h"
#include <random>
#include <map>
#include <Eigen/Core>
#include <algorithm> // std::unique_copy
#include <iostream>// writing on a text file
#include <fstream>
#include <math.h>
#include <assert.h>

using namespace std;
using namespace Aboria;
using namespace Eigen; // objects VectorXf, MatrixXf

double Kd(double d){

    double Kdist;

    Kdist = 1.0-d;//-d;
//    if (d<0.5){
//        Kdist = 1.0;
//    }else{
//        Kdist = 0.0;
//    }

    return Kdist;
}

int main() {

    int Nsim = 100;
    srand(time(0));

    for (int sim=0; sim< Nsim; sim++){


    int length_x = 1;
    int length_y = 1;
    //int n_seed = 20;

    int N = 1000; // initial number of individuals
    double beta = 0.0003;//2;//03;0.02;//
    double mu = 0.1;




    double t = 0;


    //ABORIA_VARIABLE(velocity,vdouble2,"velocity")
    ABORIA_VARIABLE(type, int, "type") // 0 - susceptible, 1 - infected, 2 - recovered
    ABORIA_VARIABLE(infpress, double, "infpress") // 0 - susceptible, 1 - infected, 2 - recovered

    typedef Particles<std::tuple<type,infpress>, 2> particle_type; // 2 stands for dimension

    // will use stored value of the position of an individual
    typedef particle_type::position position;

    // initialise the number of particles
    particle_type particles(N);

    /*
     *  initialisation of individuals
     */

    // initialise random number generator for initial positions
    std::default_random_engine genx;
    std::uniform_real_distribution<double> uniformx(0, length_x);
    std::default_random_engine geny;
    geny.seed(sim); // choose different seeds to obtain different random numbers
    std::uniform_real_distribution<double> uniformy(0, length_y);
    //srand((unsigned)time(NULL));// seed based on simulation number


    for (int i = 0; i < N; ++i) {

        double a,b;
            a = (double) rand()/RAND_MAX;
            b = (double) rand()/RAND_MAX;


        get<type>(particles[i]) = 0; // all susceptible
        get<infpress>(particles[i]) = 0; // no infectious pressure

        get<position>(particles[i]) = vdouble2(a,b);

    }

    // one particle infected
    get<type>(particles)[0] = 1; //assume that the first one since all the coordinates are random



    // initialise neighbourhood search, note that the domain will grow in x direction, so I initialise larger domain
    particles.init_neighbour_search(vdouble2(0, 0), vdouble2(length_x, length_y), vbool2(false, false));

    // save particles before they move

    //vtkWriteGrid("Initial", t, particles.get_grid(true));

    // initialise random number generator to obtain random number between 0 and 2*pi
    std::default_random_engine gen1;
    gen1.seed(sim); // choose different seeds to obtain different random numbers
    std::uniform_real_distribution<double> uniformstd(0, 1);
    std::uniform_real_distribution<double> uniformdistance(0, length_x);


    int countInf =1; // Infectious initially
    int countSus = N-countInf; // Susceptible initially
    int countRec = 0;  // Recovered initially

    //srand(time(0));

        while (countInf >0){


//        #ifdef HAVE_VTK
//                vtkWriteGrid("NewIndividualsKd1data", t*1000, particles.get_grid(true));
//        #endif



        double rand1,rand2;
        rand1 = (double) rand()/RAND_MAX;
        rand2 = (double) rand()/RAND_MAX;

        //double dist = uniformdistance(gen1);

        // for every susceptible k, calculate infectious pressure on that individual j as a sum_{infected i) beta K(d)
        VectorXd infpressure = VectorXd::Zero(countSus);
        int c=0;
        for (int k=0; k< particles.size(); k++){
            if (get<type>(particles[k]) == 0){ // for every susceptible
                for (int m=0; m< particles.size();m++){ // look for infected individuals

                    if (get<type>(particles[m]) == 1){ // if infected, add the infectious pressure
                        vdouble2 diff = get<position>(particles[k]) - get<position>(particles[m]);
                        infpressure[c] = infpressure[c] + beta*Kd(diff.norm());
                    }

                }
                get<infpress>(particles[k]) = infpressure[c];
                c = c+1;
            }

        }

        t = t + 1/ (infpressure.sum() + mu * countInf) * log(1/rand1);

//        // pick a random infectious individual
//        int rand3 = (rand() % countInf) + 1;
//        int speccountInf=0;
//        int targetInf;
//        for (int i = 0; i < particles.size(); i++) {
//            if (get<type>(particles)[i]==1) {
//                speccountInf = speccountInf + 1;
//                if (speccountInf == rand3){
//                    targetInf = i;
//                    break;
//                }
//            }
//        }
//        vdouble2 x =get<position>(particles[targetInf]);

//        if (countInf == 1){
//            cout << "rand2 " << rand2 << endl;
//            cout << "mu*countInf/(infpressure.sum()+mu*countInf)) " << mu*countInf/(infpressure.sum()+mu*countInf) << endl;
//        }


        if (rand2 < mu*countInf/(infpressure.sum()+mu*countInf)){ //recover

        // pick a random infectious individual
        int rand3 = (rand() % countInf) + 1;
        int speccountInf=0;
        int targetInf;
        for (int i = 0; i < particles.size(); i++) {
            if (get<type>(particles)[i]==1) {
                speccountInf = speccountInf + 1;
                if (speccountInf == rand3){
                    targetInf = i;
                    break;
                }
            }
        }
        get<type>(particles[targetInf]) = 2; // recovered

            // target infectious individual transmits the disease to an indivdual d distance away
//
//            for (auto k = euclidean_search(particles.get_query(), x, dist); k != false; ++k) {
//
//                if (get<type>(*k) == 0) { // if a neigbouring individual is susceptible
//                    get<type>(*k) = 1;
//                    break;
//                }
//            }
        }else{//
            // Pick a susceptible to infect based on a probability that is proportion to its infectious pressure (compared to the total infectious pressure over all individuals).

            // pick a random number
            double r;
            r = (double) rand()/RAND_MAX;
            double prob = 0;
            double prob_old = 0;
            // for each susceptible calculate it's likelihoood of being infected
            for (int p = 0;p< particles.size();p++){
                //cout << " p " << p << endl;
                if (get<type>(particles[p]) == 0){
                    //cout << "susc " << endl;
                    prob = prob + get<infpress>(particles[p])/infpressure.sum();

                    if (prob_old < r <= prob){
                        get<type>(particles[p]) = 1;
//                        cout << "prob_old " << prob_old << endl;
//                        cout << "prob" << prob << endl;
//                        cout << "changed" << endl;
                        break;

                    }
                    prob_old = prob_old + get<infpress>(particles[p])/infpressure.sum();
                }
            }
            //cout << "left " << endl;

        }


        // Count infectious
        countInf = 0;
        for (int i = 0; i < particles.size(); i++) {
                if (get<type>(particles)[i]==1){
                    countInf = countInf +1;
                }}


            // count Susceptibles
        countSus = 0;
        for (int i = 0; i < particles.size(); i++) {
            if (get<type>(particles)[i]==0){
                countSus = countSus +1;
            }
        }

    }
        // count Recovered
    for (int i = 0; i < particles.size(); i++) {
        if (get<type>(particles)[i]==2){
            countRec = countRec +1;
        }
    }

//    cout << "t " << t << endl;

    cout << countRec << endl;
    }
}




