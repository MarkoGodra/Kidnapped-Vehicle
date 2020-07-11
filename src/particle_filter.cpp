/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 150;  // TODO: Set the number of particles
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

  std::default_random_engine gen;

    // This line creates a normal (Gaussian) distribution for x
  std::normal_distribution<double> dist_x(x, std_x);
  std::normal_distribution<double> dist_y(y, std_y);
  std::normal_distribution<double> dist_theta(theta, std_theta);

  // Resize weights vector
  weights.resize(num_particles);

  for(auto i = 0; i < num_particles; i++)
  {
    // Initialize particle
    Particle particle;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;

    particles.push_back(particle);

    // Initialize coresponding weight
    weights[i] = 1.0;
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  // Unpack std deviations
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];

  // Create normal distirubtions
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(0, std_x);
  std::normal_distribution<double> dist_y(0, std_y);
  std::normal_distribution<double> dist_theta(0, std_theta);

  for(auto i = 0; i < num_particles; i++)
  {
    // Initial x,y,theta are particles current values
    double x0 = particles[i].x;
    double y0 = particles[i].y;
    double theta0 = particles[i].theta;

    // Avoid division by zero
    if(fabs(yaw_rate) > 0.001)
    {
      double v_theta_dot = velocity / yaw_rate;

      // Update with regular motion model
      particles[i].x = x0 + (v_theta_dot * (sin(theta0 + yaw_rate * delta_t) - sin(theta0)));
      particles[i].y = y0 + (v_theta_dot * (cos(theta0) - cos(theta0 + yaw_rate * delta_t)));
      particles[i].theta = theta0 + (yaw_rate * delta_t);
    }
    else
    {
      // Straight line - no change of yaw
      particles[i].x = x0 + velocity * delta_t * cos(theta0);
      particles[i].y = y0 + velocity * delta_t * sin(theta0);
    }

    // Add random noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  // TODO: For each particle
  // 1. Transform observation (car coordinate sys -> map coordinate sys)
  // 2. Associate transformed observations with nearest map landmarks
  // 3. Update weights
  // 3a. Multi dimensional gausian func
  // 3b. Combine probs

  for(size_t i = 0; i < particles.size(); i++)
  {
    // 1. Transform observations to map and then move to particle referent frame
    
    // Prepare particle coordinates
    double xp = particles[i].x;
    double yp = particles[i].y;
    double thetap = particles[i].theta;
    double weight = 1.0;

    // Prepare vector in which we shall place landmarks that are within sensor range
    std::vector<LandmarkObs> landmarks_in_range;
    
    for(const auto& landmark : map_landmarks.landmark_list)
    {
      // Only consider map landmarks that are in sensor range from current particle
      if(dist(xp, yp, landmark.x_f, landmark.y_f) < sensor_range)
      {
        LandmarkObs observed_landmark;
        observed_landmark.x = landmark.x_f;
        observed_landmark.y = landmark.y_f;
        observed_landmark.id = landmark.id_i;
        landmarks_in_range.push_back(observed_landmark);
      }
    }

    // Iterate through observations
    for(const auto& observation : observations)
    {
      double xc = observation.x;
      double yc = observation.y;

      // For each observation convert it from car coordinate system to map coordinate systems (car c.s. -> map c.s.)
      double xm = xp + cos(thetap) * xc - sin(thetap) * yc;
      double ym = yp + sin(thetap) * xc + cos(thetap) * yc;

      double min_dist_obs_to_landmark = sensor_range;
      LandmarkObs nearest_landmark;
      
      // Find landmark that is nearest to current observation

      // For each landmark
      for(const auto& landmark : landmarks_in_range)
      {
        // Measure current distance from current observation to landmark
        double dist_obs_to_landmark = dist(xm, ym, landmark.x, landmark.y);

        // If distance is smaller than previous min
        if(dist_obs_to_landmark < min_dist_obs_to_landmark)
        {
          // Update nearest landmark and comparison distance
          nearest_landmark = landmark;
          min_dist_obs_to_landmark = dist_obs_to_landmark;
        }
      }

      // Calculate gaussian distribution
      weight *= multivProb(std_landmark[0], std_landmark[1], xm, ym, nearest_landmark.x, nearest_landmark.y);
    }

    // Update particle weight
    particles[i].weight = weight;
    weights[i] = weight;
  }


}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  // Prepare vector for new set of particles
  std::vector<Particle> surviving_particles;

  std::default_random_engine gen;
  std::discrete_distribution<int> distribution(weights.begin(), weights.end());
  
  for(auto i = 0; i < num_particles; i++)
  {
    int index = distribution(gen);
    surviving_particles.push_back(particles[index]);
  }

  particles = surviving_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

double ParticleFilter::multivProb(double sig_x, double sig_y, double x_obs, double y_obs, double mu_x, double mu_y)
{
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}