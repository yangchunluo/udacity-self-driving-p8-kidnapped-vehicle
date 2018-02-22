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
#include <limits>
#include <cassert>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
  //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  default_random_engine gen;

  // Create Gaussian distributions for each positional parameter.
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  // For each particle, draw from the above distributions.
  for (int i = 0; i < this->num_particles; i++) {
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0f;
    this->particles.emplace_back(p);
  }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

  // For the random Gaussian noise, ideally it should come from std_velocity and std_yaw_rate.
  // Using std_pos is a vast simplification in this project.
  default_random_engine gen;

  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  for (auto& p : this->particles) {
    p.x += fabs(yaw_rate) < 0.000001 ?
               velocity * cos(p.theta) * delta_t :
               velocity / yaw_rate * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
    p.x += dist_x(gen);
    
    p.y += fabs(yaw_rate) < 0.000001 ?
               velocity * sin(p.theta) * delta_t :
               velocity / yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t));
    p.y += dist_y(gen);

    p.theta += yaw_rate * delta_t;
    p.theta += dist_theta(gen);
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
  //   and the following is a good resource for the actual equation to implement (look at equation 3.33)
  //   http://planning.cs.uiuc.edu/node99.html

  for (auto& p : this->particles) {

    p.weight = 1.0f;
    p.associations.clear();
    p.sense_x.clear();
    p.sense_y.clear();

    // Filter out the landmarks that is outside the sensor range for this particle.
    vector<LandmarkObs> lm_in_range;
    for (auto& l : map_landmarks.landmark_list) {
      if (dist(p.x, p.y, l.x_f, l.y_f) > sensor_range) {
        continue;
      }
      lm_in_range.emplace_back(LandmarkObs(l.id_i, l.x_f, l.y_f));
    }

    // For each observation, associate it with a landmark.
    for (auto& o : observations) {
      // Translation from vehicle (particle) coordinate system to map coordinate system.
      const double xm = p.x + cos(p.theta) * o.x - sin(p.theta) * o.y;
      const double ym = p.y + sin(p.theta) * o.x + cos(p.theta) * o.y;

      // Data association using nearest-neighbor heuristics.
      double min_dist = numeric_limits<double>::max();
      int index = -1;
      for (int i = 0; i < lm_in_range.size(); i++) {
        double curr_dist = dist(xm, ym, lm_in_range[i].x, lm_in_range[i].y);
        if (curr_dist < min_dist) {
          min_dist = curr_dist;
          index = i;
        }
      }

      // Compute the multivariate-gaussian probability for this observation.
      double prob = 0;
      if (index >= 0) {
        double lx = lm_in_range[index].x;
        double ly = lm_in_range[index].y;
        prob = exp(-(squared(xm - lx) / 2 / squared(std_landmark[0]) +
                     squared(ym - ly) / 2 / squared(std_landmark[1])))
                   / 2 / M_PI / std_landmark[0] / std_landmark[1];
      } else {
        cout<<"no landmark in range"<<endl;
      }

      // Update this particle's weight.
      p.weight *= prob;

      // Bookeeping.
      p.sense_x.emplace_back(xm);
      p.sense_y.emplace_back(ym);
      p.associations.emplace_back(index >= 0 ? lm_in_range[index].id : -1);
    } 
  }
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight. 
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  vector<double> weights;
  for (auto& p : this->particles) {
    weights.emplace_back(p.weight);
  }
  discrete_distribution<int> dist(weights.begin(), weights.end());
  default_random_engine gen;

  vector<Particle> new_particles;
  for (int i = 0; i < this->num_particles; i++) {
    new_particles.emplace_back(this->particles[dist(gen)]);
  }

  this->particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x,
										 const std::vector<double>& sense_y) {
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseX(Particle best) {
  vector<double> v = best.sense_x;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseY(Particle best) {
  vector<double> v = best.sense_y;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
