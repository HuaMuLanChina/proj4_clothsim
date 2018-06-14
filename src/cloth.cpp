#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

void Cloth::buildGrid() {
  // TODO (Part 1.1): Build a grid of masses.
	double width_step = width / num_width_points;
	double height_step = height / num_height_points;
	for(int i = 0 ; i < num_height_points; i++)
	{
		for (int j = 0; j < num_width_points; j++)
		{
			Vector3D pos(i*height_step, 0, j*width_step);
			PointMass mass(pos, false);
			point_masses.emplace_back(mass);
		}
	}
		
	for (int k = 0; k < pinned.size(); k++)
	{
		int x, y;
		x = pinned[k][0];
		y = pinned[k][1];

		point_masses[y*num_height_points+x].pinned = true;
	}
  // TODO (Part 1.2): Add springs 
	//STRUCTURAL
	for (int y = 0; y < num_height_points; y++) {
		for (int x = 0; x < num_width_points - 1; x++) {
			PointMass *pm = &point_masses[y * num_width_points + x];
			Spring structual_spring(pm, pm+1, e_spring_type::STRUCTURAL);
			springs.emplace_back(structual_spring);
		}
	}

	for (int y = 0; y < num_height_points - 1; y++)
	{
		for (int x = 0; x < num_width_points; x++)
		{
			PointMass *pm = &point_masses[y * num_width_points + x];
			Spring structual_spring(pm, pm + num_width_points, e_spring_type::STRUCTURAL);
			springs.emplace_back(structual_spring);
		}

	}
	//SHEARING
	for (int y = 0; y < num_height_points - 1; y++)
	{
		for (int x = 0; x < num_width_points - 1; x++)
		{
			PointMass *pm = &point_masses[y * num_width_points + x];
			Spring shearing_spring0(pm, pm + num_width_points + 1, e_spring_type::SHEARING);
			Spring shearing_spring1(pm + 1, pm + num_width_points, e_spring_type::SHEARING);
			springs.emplace_back(shearing_spring0);
			springs.emplace_back(shearing_spring1);
		}
	}
	//BENDING
	for (int y = 0; y < num_height_points - 2; y++)
	{
		for (int x = 0; x < num_width_points - 2; x++)
		{
			PointMass *pm = &point_masses[y * num_width_points + x];
			Spring bending_spring0(pm, pm + 2, e_spring_type::BENDING);
			Spring bending_spring1(pm, pm + 2*num_width_points, e_spring_type::BENDING);
			springs.emplace_back(bending_spring0);
			springs.emplace_back(bending_spring1);
		}
	}

}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2.1): Compute total force acting on each point mass.
  for (int y = 0; y < num_height_points; y++)
  {
	  for (int x = 0; x < num_width_points; x++)
	  {
		  PointMass *pm = &point_masses[y * num_width_points + x];
		  pm->forces = 0;
	  }
  }
  for each (Vector3D force in external_accelerations)
  {
	  for (int y = 0; y < num_height_points; y++)
	  {
		  for (int x = 0; x < num_width_points; x++)
		  {
			  PointMass *pm = &point_masses[y * num_width_points + x];
			  pm->forces += force * mass;
		  }
	  }
  }

  for each(Spring spring in springs)
  {
	  switch (spring.spring_type)
	  {
		  case e_spring_type::STRUCTURAL:
			  if (cp->enable_structural_constraints)
			  {
				  Vector3D dir = spring.pm_a->position - spring.pm_b->position;
				  double force = cp->ks*(dir.norm() - spring.rest_length);
				  Vector3D sForce = dir.unit() * force;
				  spring.pm_a->forces += -sForce;
				  spring.pm_b->forces += sForce;
			  }
			  break;
		  case e_spring_type::SHEARING:
			  if (cp->enable_structural_constraints)
			  {
				  Vector3D dir = spring.pm_a->position - spring.pm_b->position;
				  double force = cp->ks*(dir.norm() - spring.rest_length);
				  Vector3D sForce = dir.unit() * force;
				  spring.pm_a->forces += -sForce;
				  spring.pm_b->forces += sForce;
			  }
			  break;
		  case e_spring_type::BENDING:
			  if (cp->enable_structural_constraints)
			  {
				  Vector3D dir = spring.pm_a->position - spring.pm_b->position;
				  double force = 0.2*cp->ks*(dir.norm() - spring.rest_length);
				  Vector3D sForce = dir.unit() * force;
				  spring.pm_a->forces += -sForce;
				  spring.pm_b->forces += sForce;
			  }
			  break;
	  default:
		  break;
	  }
  }

  // TODO (Part 2.2): Use Verlet integration to compute new point mass positions

  for (int y = 0; y < num_height_points; y++)
  {
	  for (int x = 0; x < num_width_points; x++)
	  {
		  PointMass *pm = &point_masses[y * num_width_points + x];
		  if (!pm->pinned)
		  {
			  Vector3D acc = (pm->forces / mass);
			  Vector3D newpos = pm->position + (1 - cp->damping)*(pm->position - pm->last_position) + acc*delta_t*delta_t;
			  pm->last_position = pm->position;
			  pm->position = newpos;
		  }
	  }
  }

  // This won't do anything until you complete Part 4.
  build_spatial_map();
  for (PointMass &pm : point_masses) {
    self_collide(pm, simulation_steps);
  }

  // This won't do anything until you complete Part 3.
  for (PointMass &pm : point_masses) {
    for (CollisionObject *co : *collision_objects) {
      co->collide(pm);
    }
  }


  // TODO (Part 2.3): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
  for each(Spring spring in springs)
  {
	  double constrain_length = spring.rest_length * 1.1;

	  double current_length = (spring.pm_a->position - spring.pm_b->position).norm();
	  if (current_length >= constrain_length)
	  {
		  Vector3D dir = (spring.pm_a->position - spring.pm_b->position).unit();
		  
		  if (!(spring.pm_a->pinned || spring.pm_b->pinned))
		  {
			  double half_correction = (current_length - constrain_length)*0.5;
			  spring.pm_a->position = spring.pm_a->position - dir * half_correction;
			  spring.pm_b->position = spring.pm_b->position + dir * half_correction;
		  }
		  else if (spring.pm_a->pinned)
		  {
			  double correction = (current_length - constrain_length);
			  spring.pm_b->position = spring.pm_b->position + dir * correction;
		  }
		  else if (spring.pm_b->pinned)
		  {
			  double correction = (current_length - constrain_length);
			  spring.pm_a->position = spring.pm_a->position - dir * correction;
		  }
	  }
  }

}

void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4.2): Build a spatial map out of all of the point masses.

}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4.3): Handle self-collision for a given point mass.

}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4.1): Hash a 3D position into a unique float identifier that represents
  // membership in some uniquely identified 3D box volume.

  return 0.f;
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm, pm + num_width_points, pm + 1));
      triangles.push_back(new Triangle(pm + 1, pm + num_width_points,
                                       pm + num_width_points + 1));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}
