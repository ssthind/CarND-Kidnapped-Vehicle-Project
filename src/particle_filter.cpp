/*
 * particle_filter.cpp
 *
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Setting the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  num_particles =100;
  default_random_engine gen; 
   
 	// This line creates a normal (Gaussian) distribution for x,y and theta.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; ++i) {
      Particle particle1;
      
      particle1.x = dist_x(gen);
      particle1.y = dist_y(gen);
      particle1.theta = dist_theta(gen);
      particle1.id=i;
      particle1.weight = 1.0;
      particles.push_back(particle1);
      weights.push_back(1.0);
      
  }
  
  is_initialized = true;
  
  
  
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Added measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  default_random_engine gen;
  double predicted_x, predicted_y, predicted_theta;
  for (unsigned int i =0; i< num_particles; i++){
    if (yaw_rate==0){
      predicted_x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
      predicted_y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
      predicted_theta = particles[i].theta;
    }
    else{
      predicted_x = particles[i].x + (velocity / yaw_rate)*( sin(particles[i].theta + delta_t * yaw_rate) - sin(particles[i].theta) );
      predicted_y = particles[i].y + (velocity / yaw_rate)*( cos(particles[i].theta) - cos(particles[i].theta + delta_t * yaw_rate) );
      predicted_theta = particles[i].theta + (delta_t * yaw_rate);      

    }
	normal_distribution<double> dist_x(predicted_x, std_pos[0]);
	normal_distribution<double> dist_y(predicted_y, std_pos[1]);
	normal_distribution<double> dist_theta(predicted_theta, std_pos[2]);
	particles[i].x = dist_x(gen);
	particles[i].y = dist_y(gen);
	particles[i].theta = dist_theta(gen);

  }
  
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// Finding the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  double prev_min_dist, next_dist; 
  for (unsigned int i=0; i<observations.size() ; i++){
	prev_min_dist = std::numeric_limits<double>::max();
    for (unsigned int j=0; j<predicted.size() ; j++){
      next_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      if (prev_min_dist>=next_dist){
        prev_min_dist=next_dist;
        observations[i].id=j;
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// Updating the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
  

  double theta_par, x_par, y_par; 
  vector<LandmarkObs> LandMark_inrange;
  LandmarkObs  obs_map, land_map;
  vector<LandmarkObs> obs_landmark_maped;
  double sum_weights;
  for(unsigned int i =0; i< num_particles; i++){
    x_par = particles[i].x;
    y_par = particles[i].y;
    theta_par = particles[i].theta;    
    sum_weights= 0.0;
    // finding Land Marks in range of particle i in consideration from list of map_landmarks
    LandMark_inrange.clear();
    for(unsigned int j =0; j< map_landmarks.landmark_list.size(); j++){
		land_map.id = map_landmarks.landmark_list[j].id_i;
		land_map.x = map_landmarks.landmark_list[j].x_f;
		land_map.y = map_landmarks.landmark_list[j].y_f;
      if (dist(x_par, y_par, land_map.x, land_map.y) <= sensor_range){
        LandMark_inrange.push_back(land_map);  
      }
    }
    //Converting all observations to map coordinates for each particle i in consideration
    obs_landmark_maped.clear();
    for(unsigned int o=0; o<observations.size() ; o++){
      obs_map.id = observations[o].id;
      obs_map.x = x_par + observations[o].x * cos(theta_par) - observations[o].y * sin(theta_par);
      obs_map.y = y_par + observations[o].x * sin(theta_par) + observations[o].y * cos(theta_par);     
      obs_landmark_maped.push_back(obs_map);      
    }
    
    // Calling landmark association function to get landmark data
    dataAssociation(LandMark_inrange, obs_landmark_maped);
    
    //associating all observerd landmarks to map
    vector<int> association1;
    vector<double> sense_x1;
    vector<double> sense_y1;
    for (unsigned int o=0; o<obs_landmark_maped.size();o++){
      association1.push_back(LandMark_inrange[obs_landmark_maped[o].id].id);
	  sense_x1.push_back(obs_landmark_maped[o].x);
      sense_y1.push_back(obs_landmark_maped[o].y);
    }
	//Associating weights
    SetAssociations(particles[i], association1, sense_x1, sense_y1);
    
    // calculate Multivariate-Gaussian probability density for weight assigment
    double weight_exp =0.0;
    double counter_exp=obs_landmark_maped.size();
    double P_constant= 1.0/(2.0* M_PI * std_landmark[0] *std_landmark[1]);
    double exp_i, delta_x, delta_y;
    for(unsigned int o=0; o<counter_exp ; o++){
      delta_x =  LandMark_inrange[obs_landmark_maped[o].id].x - obs_landmark_maped[o].x;
      delta_y =  LandMark_inrange[obs_landmark_maped[o].id].y - obs_landmark_maped[o].y;
      exp_i = -((delta_x*delta_x )/(2.0*std_landmark[0] * std_landmark[0]) + (delta_y*delta_y )/(2.0*std_landmark[1] * std_landmark[1]));
	  //Check if exponent value is getting to high, to avoid zero probablity value
      if (fabs(exp_i)<100){
      weight_exp += exp_i;  
      }
      else{
      weight_exp = -100;
      }
    }
	
	//weights calculation
	weights[i]= pow(P_constant, counter_exp) * exp(weight_exp);
    particles[i].weight = weights[i];
	// finding sum of weights for normalization
    sum_weights+=particles[i].weight;
	}
	// normalization weights
	for(unsigned int i =0; i< num_particles; i++){
		weights[i] /= sum_weights;
		particles[i].weight = weights[i];
	}
}

void ParticleFilter::resample() {
	// Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<int> idx(weights.begin(), weights.end());
  std::vector<Particle> new_particles;
  
  for(unsigned int i =0; i< num_particles; i++){
    new_particles.push_back(particles[idx(gen)]);
  }
  
  particles = new_particles;

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
	
	return particle;
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
