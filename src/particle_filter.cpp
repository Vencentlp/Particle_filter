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
    num_particles = 50;
   // is_initialized = false;
    
    //Define ranom engine
    default_random_engine gen;
    
    
    //Creat normal distribution for x,y,theta
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    
    //sample from the distribution above to creat particles
    
    Particle single_particle;
    
    for (int i=0; i<num_particles; i++)
    {
        single_particle.x = dist_x(gen);
        single_particle.y = dist_y(gen);
        single_particle.theta = dist_theta(gen);
        particles.push_back(single_particle);
        //std::cout << single_particle.x << std::endl;
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
        
        particles[i].x = x + dist_x(gen);
        particles[i].y = y + dist_y(gen);
        while(theta < -M_PI){theta += 2*M_PI;}
        while(theta > M_PI){theta -=2*M_PI;}
        particles[i].theta = theta + dist_theta(gen);
        
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
            double ym = xc * sin(theta) + yc * cos(theta) + yc;
            single_translandmark.x = xm;
            single_translandmark.y = ym;
            single_translandmark.id = observations[j].id;
            trans_landmarks.push_back(single_translandmark);                       
        }
        
        for (int k=0;k<map_landmarks.landmark_list.size();k++)
        {
            LandmarkObs single_landmark;
            double xf = map_landmarks.landmark_list[k].x_f;
            double yf = map_landmarks.landmark_list[k].y_f;
            double dist_r = dist(xf,yf,xp,yp);
            
            
            
            if (dist_r<=sensor_range)
            {
                single_landmark.x = map_landmarks.landmark_list[k].x_f;
                single_landmark.y = map_landmarks.landmark_list[k].y_f;
                single_landmark.id = map_landmarks.landmark_list[k].id_i;
                predicted.push_back(single_landmark);
            }
        }
        
        
        dataAssociation(predicted, trans_landmarks);
        
        double weight = 1;
        for (int i=0;i<trans_landmarks.size();i++)
        {
            int id = trans_landmarks[i].id;
            double expox = trans_landmarks[i].x-map_landmarks.landmark_list[id].x_f;
            double expoy = trans_landmarks[i].x-map_landmarks.landmark_list[id].x_f;
            double expo = -expox*expox/(2*std_landmark[0]*std_landmark[0])-expoy*expoy/(2*std_landmark[1]*std_landmark[1]);
            weight *= 0.5/(M_PI*std_landmark[0]*std_landmark[1])*exp(expo);
        }
        weights.push_back(weight);
        
    }
    
    //Sum of the weights
    double sum = accumulate(weights.begin(), weights.end(),0);
    //Normalize particle weights
    for (int i=0;i<particles.size();i++)
    {
        particles[i].weight = weights[i]/sum;
        weights[i] = weights[i]/sum;
    } 

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    double beta = 0;
    vector<Particle> particles_resample;
    default_random_engine gen;
    uniform_int_distribution<int> index_dist(0,num_particles-1);
    auto weight_max = max_element(weights.begin(),weights.end());
    uniform_real_distribution<double> beta_dist(0,2 * (*weight_max));
    
    int index = index_dist(gen);
    
    for (int i=0;i<num_particles;i++)
    {
        beta = beta + beta_dist(gen);
        while(weights[index] < beta)
        {
            beta = beta -weights[index];
            index = (index + 1)%num_particles;            
        }
        particles_resample.push_back(particles[index]);
        
    }
    particles = particles_resample;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
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
