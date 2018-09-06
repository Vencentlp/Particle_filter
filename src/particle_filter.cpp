/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>      
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"
//#include "helper_functions.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles = 100;
   // is_initialized = false;
    
    //Define ranom engine
    default_random_engine gen;
    
    weights.resize(num_particles);
    
    
    //Creat normal distribution for x,y,theta
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    //std::cout<<"GPS"<<x<<" "<<y<<" "<<theta<<std::endl;
    
    //sample from the distribution above to creat particles
    
    Particle single_particle;
    
    for (int i=0; i<num_particles; i++)
    {
        single_particle.x = dist_x(gen);
        single_particle.y = dist_y(gen);
        single_particle.theta = dist_theta(gen);
        particles.push_back(single_particle);
        //std::cout << "Particle:"<< i << std::endl<<single_particle.x<<" "<<single_particle.y<<" "<<single_particle.theta<< std::endl;
    }
    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    //define the random engine
    default_random_engine gen;
    
    //Predict for each particle
    for (int i=0;i<particles.size();i++)
    {
        double theta_prev = particles[i].theta;
        double x_prev = particles[i].x;
        double y_prev = particles[i].y;
        double x,y,theta;
        //std::cout << "Particle_previous:"<< i << std::endl<<x_prev<<" "<<y_prev<<" "<<theta_prev<< std::endl;
        //std::cout<<"deltat:"<<delta_t<<std::endl;
        
        
        if (abs(yaw_rate)<0.0001)
        {
            x = x_prev + velocity*cos(theta_prev)*delta_t;
            y = y_prev + velocity*sin(theta_prev)*delta_t;            
        }
        else
        {
            x = x_prev + velocity/yaw_rate*(sin(theta_prev+yaw_rate*delta_t) - sin(theta_prev));
            y = y_prev + velocity/yaw_rate*(-cos(theta_prev+yaw_rate*delta_t)+cos(theta_prev));             
        }
        theta = theta_prev + yaw_rate*delta_t;
        //Creat normal distribution for x,y,theta
        normal_distribution<double> dist_x(x, std_pos[0]);
        normal_distribution<double> dist_y(y, std_pos[1]);
        normal_distribution<double> dist_theta(theta, std_pos[2]);        

        //std::cout << "Dist:"<< i << std::endl<<x_dist<<" "<<y_dist<<" "<<theta_dist<< std::endl;
        particles[i].x = dist_x(gen);;
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
        //while(particles[i].theta < -M_PI){particles[i].theta += 2*M_PI;}
        //while(particles[i].theta > M_PI){particles[i].theta -=2*M_PI;}
         
        //std::cout << "Motion control"<<velocity<<" "<<theta<<std::endl;
        //std::cout<<"observations:"<<observations[j].id<<" "<<observations[j].x<<" "<<observations[j].y<<std::endl;
        
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    
    for (int i=0;i<observations.size();i++)
    {
        double sub = dist(observations[i].x,observations[i].y,predicted[0].x,predicted[0].y);
        for (int j=0;j<predicted.size();j++)
        {
            double dist_l = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
            if (sub >= dist_l)
            {
                sub = dist_l;
                observations[i].id = predicted[j].id;                
            }
            
        }
       // std::cout<<"dist"<<sub<<endl;
    }
    
    

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html   
    
    //Transform the observations coordinate from vehicle coordinates into map coordinates
    for (int i=0;i<particles.size();i++)
    {
        LandmarkObs single_translandmark;
        vector<LandmarkObs> trans_landmarks;
        vector<LandmarkObs> predicted;
        
        
        double xp = particles[i].x;
        double yp = particles[i].y;
        double theta = particles[i].theta;
        
        //transform the observation into map coordinate
        for (int j=0; j<observations.size();j++)
        {
            double xc = observations[j].x;
            double yc = observations[j].y;
            double xm = xc * cos(theta) - yc * sin(theta) + xp;
            double ym = xc * sin(theta) + yc * cos(theta) + yp;
            single_translandmark.x = xm;
            single_translandmark.y = ym;
            single_translandmark.id = observations[j].id;
            trans_landmarks.push_back(single_translandmark); 
            
           // std::cout << "Particle:"<< i << std::endl<<xp<<" "<<yp<<" "<<theta<< std::endl;
           // std::cout<<"observations:"<<observations[j].id<<" "<<observations[j].x<<" "<<observations[j].y<<std::endl;
           // std::cout<<"trans_landmarks:"<<trans_landmarks[j].id<<" "<<trans_landmarks[j].x<<" "<<trans_landmarks[j].y<<std::endl;
        }
        
        //Extract the landmarks with the sensor range and stroed in predicted vector
        for (int k=0;k<map_landmarks.landmark_list.size();k++)
        {
            LandmarkObs single_landmark;
            double xf = map_landmarks.landmark_list[k].x_f;
            double yf = map_landmarks.landmark_list[k].y_f;
            double dist_r = dist(xf,yf,xp,yp);
            
            
            //double thresh = 5.0;
            if (dist_r<=(sensor_range))
            {
                single_landmark.x = map_landmarks.landmark_list[k].x_f;
                single_landmark.y = map_landmarks.landmark_list[k].y_f;
                single_landmark.id = map_landmarks.landmark_list[k].id_i;
                predicted.push_back(single_landmark);           
                
            }
        }
        
        
        
        // Associate the closest landmarks with each observation
        dataAssociation(predicted, trans_landmarks);
       
        /*+for (int p=0;p<predicted.size();p++)
        {
            std::cout<<"Predicted_landmarks:"<<predicted[p].id<<" "<<predicted[p].x<<" "<<predicted[p].y<<std::endl;
        }
        
        for (int p=0;p<trans_landmarks.size();p++)
        {
            std::cout<<"Trans_landmarks:"<<trans_landmarks[p].id<<" "<<trans_landmarks[p].x<<" "<<trans_landmarks[p].y<<std::endl;
        }
        
         for (int p=0;p<map_landmarks.landmark_list.size();p++)
        {
            std::cout<<"map_landmarks:"<<map_landmarks.landmark_list[p].id_i<<" "<<map_landmarks.landmark_list[p].x_f<<" "<<map_landmarks.landmark_list[p].y_f<<std::endl;
        }
        */
        
        
        //Calculate the multigaussian weight
        double weight = 1;
        vector<int> associations;
        vector<double> sense_x;
        vector<double> sense_y;
        for (int l=0;l<trans_landmarks.size();l++)
        {
            int id = trans_landmarks[l].id;
            double x = trans_landmarks[l].x;
            double y = trans_landmarks[l].y;
            sense_x.push_back(x);
            sense_y.push_back(y);
            associations.push_back(id);
            
            double expox = x-map_landmarks.landmark_list[id-1].x_f;
            double expoy = y-map_landmarks.landmark_list[id-1].y_f;
            double expo = -expox*expox/(2*std_landmark[0]*std_landmark[0])-expoy*expoy/(2*std_landmark[1]*std_landmark[1]);
            //double p_gaus = 0.5/(M_PI*std_landmark[0]*std_landmark[1])*exp(expo);
            //if (abs(p_gaus) > 0.0001)
            //{
            weight *= 0.5/(M_PI*std_landmark[0]*std_landmark[1])*exp(expo);
            //}
            //std::cout<<"trans_landmarks_afterass:"<<trans_landmarks[l].id<<" "<<trans_landmarks[l].x<<" "<<trans_landmarks[l].y<<std::endl;
            //std::cout<<"num"<<l<<" expox: "<<expox<<" expoy"<<expoy<<" expo"<<expo<<" weight"<<weight<<endl;
            
        }
        weights[i] = weight;
        
        //Set associations
        SetAssociations(particles[i],associations,sense_x,sense_y);
        
        //std::cout<<"particle"<<i<<"weight"<<weight<<"weightssize"<<weights.size()<<endl;
        
    }
    
    //Sum of the weights
    double sum = 0;
    for (unsigned int i=0;i<weights.size();i++){sum+=weights[i];}
    std::cout<<"sum"<<": "<<sum<<endl;
    
    //Normalize particle weights
    for (unsigned int m=0;m<particles.size();m++)
    {
        particles[m].weight = weights[m]/sum;
        weights[m] = weights[m]/sum;
        //std::cout<<"Particle"<<m<<std::endl<<"weight"<< weights[m]<<std::endl;
    } 

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    double beta = 0;
    //vector<Particle> particles_resample;
    default_random_engine gen;
    //Creat a distribution for index
    uniform_int_distribution<int> index_dist(0,num_particles-1);
    //Get the max value of the weights
    auto weight_max = max_element(weights.begin(),weights.end());
    
    //std::cout<<"weightmax"<<(*weight_max)<<endl;
    
    //Creat uniform distribution for beta
    uniform_real_distribution<double> beta_dist(0,2 * (*weight_max));
    
    //Random init value fot index
    int index = index_dist(gen);
    
   //Resample using the wheel method
    for (int i=0;i<num_particles;i++)
    {
        beta = beta_dist(gen);
        //std::cout<<"Beta"<<beta<<endl;
        while(weights[index] < beta)
        {
            beta = beta -weights[index];
            index = (index + 1)%num_particles;            
        }
        particles[i] = particles[index];
        
    }
     

}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
