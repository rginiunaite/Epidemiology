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
    double a = 1.7; // 20m since 2km is 2 at this scale

    Kdist = exp(-d/a) * 1/(2*M_PI*a*a);

    //Kdist = 1.0;


    return Kdist;
}

int main() {

    int Nsim = 1;
    srand(time(0));

    bool infect_infected = false; // if infected can be infected
    bool orchards = true; // true if orcihds, false if random
    VectorXd betavalues = VectorXd::Zero(Nsim);//if only susceptibles


    //#pragma omp parallel for
    for (int sim=0; sim< Nsim; sim++){



    //int n_seed = 20;

    int N=2000;

    //double beta = 0.0003;//2;//03;0.02;// it won't be used, I will calculate beta
    double mu = 0.1;

//    int length_x = 2;
//    int length_y = 2;

    // for orchards
    double rowspacing = 10.0;
    double columnspacing = 5.0;
    int columndim = 63;
    int rowdim = 32;
    int length_x = int(rowspacing)*rowdim;
    int length_y = int(columndim)*columndim;
    if (orchards == true){
       N = columndim*rowdim; // initial number of individuals
    }

    double t = 0;


    //ABORIA_VARIABLE(velocity,vdouble2,"velocity")
    ABORIA_VARIABLE(type, int, "type") // 0 - susceptible, 1 - infected, 2 - recovered
    ABORIA_VARIABLE(infpress, double, "infpress") // 0 - susceptible, 1 - infected, 2 - recovered

    typedef Particles<std::tuple<type,infpress>, 2> particle_type; // 2 stands for dimension

    // will use stored value of the position of an individual
    typedef particle_type::position position;

    // initialise the number of particles
    particle_type particles(N);

    if (orchards == true){
        particle_type particles(columndim*rowdim);
    }


        /*
         *  initialisation of individuals
         */



    // orchards
    if (orchards == true){
        double a,b;
        int n=0;
        for (int i = 0; i < columndim; ++i) {
            for(int j = 0; j<rowdim; ++j){

                a = columnspacing* double(i);
                b = rowspacing*double(j);
                get<type>(particles[i]) = 0; // all susceptible
                get<infpress>(particles[i]) = 0; // no infectious pressure

                get<position>(particles[n]) = vdouble2(a,b);
                n=n+1;
            }

        }
    }else{
        // random initialisation
    for (int i = 0; i < N; ++i) {

        double a,b;
            a = (double) rand()/RAND_MAX;
            b = (double) rand()/RAND_MAX;


        get<type>(particles[i]) = 0; // all susceptible
        get<infpress>(particles[i]) = 0; // no infectious pressure

        get<position>(particles[i]) = vdouble2(a,b);

    }
    }



    // one particle infected
    if (orchards == true){
        double arand;
        arand = (double) rand()/RAND_MAX;
        int randnr = floor (arand * double(N)); //choose a random individual
        get<type>(particles)[randnr] = 1; //assume that the first one since all the coordinates are random
        cout << " infected position " << get<position>(particles)[randnr] << endl;

    }else{
        get<type>(particles)[0] = 1; //assume that the first one since all the coordinates are random
    }




    // initialise neighbourhood search, note that the domain will grow in x direction, so I initialise larger domain
    particles.init_neighbour_search(vdouble2(0, 0), vdouble2(length_x, length_y), vbool2(false, false));

    // save particles before they move

    //vtkWriteGrid("Initial", t, particles.get_grid(true));

    // initialise random number generator to obtain random number between 0 and 2*pi
    std::default_random_engine gen1;
    gen1.seed(sim); // choose different seeds to obtain different random numbers
    std::uniform_real_distribution<double> uniformstd(0, 1);
    //std::uniform_real_distribution<double> uniformdistance(0, length_x);


    int countInf = 1; // Infectious initially
    int countSus = N-countInf; // Susceptible initially
    int countRec = 0;  // Recovered initially



    //find beta particular to some K(d)
    VectorXd infpressureIn = VectorXd::Zero(N);//if only susceptibles
    int c= 0;
    int R0=3;

    for (int k=0; k< particles.size(); k++){
        for (int m=0; m< particles.size();m++){ // look for infected individuals
            if (get<id>(particles[k]) != get<id>(particles[m])){
                vdouble2 diff = get<position>(particles[k]) - get<position>(particles[m]);
                infpressureIn[c] = infpressureIn[c] + Kd(diff.norm());
            }
        }
        c=c+1;
    }

    double betanew = R0 * mu/infpressureIn.mean();

    betavalues[sim] = betanew;

    //cout << "betanew " << betanew << endl;




    //srand(time(0));

        //while (countInf > 0){
        while (countInf > 0 && countInf+countRec < 100){

        #ifdef HAVE_VTK
                vtkWriteGrid("O5rchids", t*1000, particles.get_grid(true));
        #endif



        double rand1,rand2;
        rand1 = (double) rand()/RAND_MAX;
        rand2 = (double) rand()/RAND_MAX;

        //double dist = uniformdistance(gen1);

        // for every susceptible k, calculate infectious pressure on that individual j as a sum_{infected i) beta K(d)
        //VectorXd infpressure = VectorXd::Zero(countSus);//if only susceptibles
        VectorXd infpressure = VectorXd::Zero(N);//if only susceptibles

            int c=0;
        for (int k=0; k< particles.size(); k++){
            //uncomment if only for susceptibles

            if (infect_infected == true){
                //if (get<type>(particles[k]) == 0 || get<type>(particles[k]) == 1){ // for every susceptible or infected
                    for (int m=0; m< particles.size();m++){ // look for infected individuals

                        if (get<type>(particles[m]) == 1){ // if infected, add the infectious pressure
                            // check if it is not the same individual
                            if (get<id>(particles[k]) != get<id>(particles[m])){
                                vdouble2 diff = get<position>(particles[k]) - get<position>(particles[m]);
                                infpressure[c] = infpressure[c] + betanew*Kd(diff.norm());
                            }
                        }

                    }
                    get<infpress>(particles[k]) = infpressure[c];
                    //                if ( get<type>(particles[k]) == 1){ // for checking
                    //                    cout << "inf pressure of infected " << infpressure[c] << endl;
                    //                }
                    c = c+1;
                    //uncomment if only for susceptibles
                //}
            }else{
                    if (get<type>(particles[k]) == 0 ){
                        for (int m=0; m< particles.size();m++){ // look for infected individuals

                            if (get<type>(particles[m]) == 1){ // if infected, add the infectious pressure
                                    vdouble2 diff = get<position>(particles[k]) - get<position>(particles[m]);
                                    infpressure[c] = infpressure[c] + betanew*Kd(diff.norm());

                            }

                        }
                        get<infpress>(particles[k]) = infpressure[c];
                        cout << "type " << get<type>(particles[k]) << endl;
                        cout << "susc postition " << get<position>(particles[k]) << endl;
                        cout << "susc infpress " << get<infpress>(particles[k]) << endl;


                        //                if ( get<type>(particles[k]) == 1){ // for checking
                        //                    cout << "inf pressure of infected " << infpressure[c] << endl;
                        //                }
                        c = c+1;
                    }
                }

        }

        cout << "infpressure.sum() " << infpressure.sum() << endl;


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
            double prob = 0.0;
            double prob_old = 0.0;
            // for each susceptible calculate it's likelihoood of being infected
            for (int p = 0;p< particles.size();p++){
                //cout << " p " << p << endl;

                if (infect_infected == true){
                    // everyone can be reinfected
                    //if (get<type>(particles[p]) == 1 || get<type>(particles[p]) == 0){

                        //cout << "susc " << endl;
                        prob = prob + get<infpress>(particles[p])/infpressure.sum();

                        if (prob_old < r <= prob){

                            // if recovered stay recovered
                            if (get<type>(particles[p]) == 2){
                                get<type>(particles[p]) = 2;
                            }else{
                                get<type>(particles[p]) = 1;
                            }
                            //                        cout << "prob_old " << prob_old << endl;
                            //                        cout << "prob" << prob << endl;
                            //                        cout << "changed" << endl;
                            break;

                        }
                        prob_old = prob_old + get<infpress>(particles[p])/infpressure.sum();

                        // uncomment the line below if only susceptibles can be infected:
                    //}
                }else{
                    // uncomment the line below if only susceptibles can be infected:
                    if (get<type>(particles[p]) == 0){

                        //cout << "susc " << endl;
                        prob = prob + double(get<infpress>(particles[p])/infpressure.sum());
                        cout << "posi " << get<position>(particles[p]) << endl;
                        cout << "prob_old " << prob_old << endl;
                        cout << "prob " << prob << endl;
                        cout << "r " << r << endl;

                        if (prob_old < r && r <= prob){
                            cout << "enters " << endl;

                            get<type>(particles[p]) = 1;
                            //                        cout << "prob_old " << prob_old << endl;
                            //                        cout << "prob" << prob << endl;
                            //                        cout << "changed" << endl;

                            cout << "infected sup position " << get<position>(particles[p]) << endl;
                            break;
                            cout << "after break " << endl;
                        }
                        prob_old = prob_old + get<infpress>(particles[p])/infpressure.sum();

                        // uncomment the line below if only susceptibles can be infected:
                    }
                }




            }
            //cout << "left " << endl;

        }


        // Count infectious
        countInf = 0;
        countSus = 0;
        countRec = 0;
        for (int i = 0; i < particles.size(); i++) {
                if (get<type>(particles)[i]==1){
                    countInf = countInf +1;
                }
            if (get<type>(particles)[i]==0){
                countSus = countSus +1;
            }
            if (get<type>(particles)[i]==2){
                countRec = countRec +1;
            }
        }

//            // count Recoverd Temporal
//            countRec = 0;
//            for (int i = 0; i < particles.size(); i++) {
//                if (get<type>(particles)[i]==2){
//                    countRec = countRec +1;
//                }
//            }
//            cout << "Recovered " <<  countRec << endl;
//            cout << "Infeced " <<  countInf << endl;
//            cout << "Susc " <<  countSus << endl;

    }
//        // count Recovered
//    for (int i = 0; i < particles.size(); i++) {
//        if (get<type>(particles)[i]==2){
//            countRec = countRec +1;
//        }
//    }

//    cout << "t " << t << endl;

    cout << countRec << endl;
    }
//
//    ofstream output("Betavalues100plantsKdexp0p1InfInfectedFalse.csv");
//
//    output << betavalues << endl;

}




