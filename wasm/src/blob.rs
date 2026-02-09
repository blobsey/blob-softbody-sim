// Blob constants
const BLOB_SPRING_STIFFNESS: f32 = 960.0;
const BLOB_SHAPE_STIFFNESS: f32 = 480.0;
const BLOB_BOUNCINESS: f32 = 0.1; // Sane values are 0.0 - 1.0
const BLOB_RADIUS: f32 = 100.0;
const BLOB_PARTICLE_RADIUS: f32 = 12.0;
const BLOB_MASS: f32 = 32.0;
const BLOB_OUTLINE_THICKNESS: f32 = 24.0;
const BLOB_SMOOTHING_PASSES: usize = 2;
const GRAVITY: f32 = 3600.0;
const BLOB_SPRING_DAMPING: f32 = BLOB_SPRING_STIFFNESS / 1.5;
const VELOCITY_DAMPING: f32 = 0.99;
const BLOB_MAX_SPEED: f32 = 100.0;
use macroquad::{
    color::{BLACK, BLUE, GREEN, RED},
    math::Vec2,
    shapes::{draw_circle, draw_line},
    window::{screen_height, screen_width},
};
use std::f32::consts::PI;

const DEBUG: bool = false;

// Blob uses a spring system for soft-body physics
// https://en.wikipedia.org/wiki/Spring_system
pub struct Blob {
    particles: Vec<Particle>,
    springs: Vec<Spring>,
    outline_particles_indices: Vec<usize>, // For drawing the outline
    center_particle_index: usize, // Nice to have as a quick center-of-mass
}

struct Particle {
    pos: Vec2,
    prev_pos: Vec2, // For Verlet integration
    // The particles "home", i.e. distance from the center of mass
    offset_from_center: Vec2,
}

struct Spring {
    particle_a: usize,
    particle_b: usize,
    rest_length: f32,
}

impl Blob {
    pub fn new(origin: Vec2) -> Blob {
        // Space between particles is 2 * radius, plus a little buffer
        let particle_spacing = BLOB_PARTICLE_RADIUS * 4.0;
        let num_outline =
            (2.0 * PI * BLOB_RADIUS / particle_spacing).floor() as usize;

        let mut positions: Vec<Vec2> = Vec::new();

        // Center particle (index 0)
        positions.push(Vec2::ZERO);
        let center_particle_index = 0;

        // Outline particles (indices 1..=num_outline)
        for i in 0..num_outline {
            let angle = 2.0 * PI * i as f32 / num_outline as f32;
            positions.push(Vec2::new(
                BLOB_RADIUS * angle.cos(),
                BLOB_RADIUS * angle.sin(),
            ));
        }

        let outline_particle_indices: Vec<usize> = (1..=num_outline).collect();

        let particles: Vec<Particle> = positions
            .iter()
            .map(|&pos| Particle {
                pos: origin + pos,
                prev_pos: origin + pos,
                offset_from_center: pos,
            })
            .collect();

        let mut springs: Vec<Spring> = Vec::new();

        // Spokes: center -> each outline particle
        for i in 0..num_outline {
            let outline_idx = 1 + i;
            springs.push(Spring {
                particle_a: center_particle_index,
                particle_b: outline_idx,
                rest_length: (particles[0].pos - particles[outline_idx].pos)
                    .length(),
            });
        }

        // Rim: adjacent outline particles
        for i in 0..num_outline {
            let curr = 1 + i;
            let next = 1 + (i + 1) % num_outline;
            springs.push(Spring {
                particle_a: curr,
                particle_b: next,
                rest_length: (particles[curr].pos - particles[next].pos)
                    .length(),
            });
        }

        Blob {
            particles,
            springs,
            outline_particles_indices: outline_particle_indices,
            center_particle_index,
        }
    }

    pub fn update(&mut self, dt: f32) {
        let mut forces = vec![Vec2::ZERO; self.particles.len()];

        // Spring forces
        for spring in &self.springs {
            let particle_a = &self.particles[spring.particle_a];
            let particle_b = &self.particles[spring.particle_b];

            let spring_vec = particle_a.pos - particle_b.pos;
            let spring_len = spring_vec.length();
            if spring_len > 0.0 {
                let unit_vec = spring_vec / spring_len;

                // Hooke's law, Force = stiffness * displacement
                let displacement = spring_len - spring.rest_length;
                let force = BLOB_SPRING_STIFFNESS * displacement;
                let force_vec = unit_vec * force;

                // Damp spring forces according to relative motion
                // i.e. if the particles are moving far apart or closer together
                // applies damping but doesn't damp if the whole blob moving
                let velocity_a = particle_a.pos - particle_a.prev_pos;
                let velocity_b = particle_b.pos - particle_b.prev_pos;
                let relative_velocity = velocity_b - velocity_a;

                // Project the velocity onto the spring axis via the dot product
                // i.e. only damp motion that compresses/stretches the string, not
                // the particles "sliding" perpendicularly
                let spring_length_change_rate = relative_velocity.dot(unit_vec);
                let damping_vec = unit_vec
                    * (BLOB_SPRING_DAMPING * spring_length_change_rate);

                forces[spring.particle_a] -= force_vec - damping_vec;
                forces[spring.particle_b] += force_vec - damping_vec;
            }
        }

        // Move each particle back towards its original position relative to
        // the "center", where center is all the particle positions averaged
        let center = self
            .particles
            .iter()
            .fold(Vec2::ZERO, |acc, particle| acc + particle.pos)
            / self.particles.len() as f32;

        for (i, particle) in self.particles.iter().enumerate() {
            let target = center + particle.offset_from_center;
            let displacement = target - particle.pos;
            forces[i] += displacement * BLOB_SHAPE_STIFFNESS;
        }

        // Gravity
        let particle_mass = BLOB_MASS / self.particles.len() as f32;

        for i in 0..forces.len() {
            forces[i].y += GRAVITY * particle_mass;
        }

        // Apply all forces
        for (i, particle) in self.particles.iter_mut().enumerate() {
            // acceleration = F/m, needed for Verlet integration
            let acceleration = forces[i] / particle_mass;

            // Verlet integration: Pₙ₊₁ = 2Pₙ - Pₙ₋₁ + accel * dt²
            let next_pos =
                2.0 * particle.pos - particle.prev_pos + acceleration * dt * dt;

            // Apply velocity damping to simulate friction
            let velocity = next_pos - particle.pos;
            let damped_velocity = velocity * VELOCITY_DAMPING;
            let damped_next_pos = particle.pos + damped_velocity;

            particle.prev_pos = particle.pos;
            particle.pos = damped_next_pos;

            let velocity = particle.pos - particle.prev_pos;
            let speed = velocity.length();
            if speed > 0.0 {
                // S-curve using tanh for smooth limiting
                let normalized_speed = speed / BLOB_MAX_SPEED;
                let s_curve_factor = normalized_speed.tanh();
                let limited_speed = s_curve_factor * BLOB_MAX_SPEED;

                let limited_velocity = velocity.normalize() * limited_speed;
                particle.prev_pos = particle.pos - limited_velocity;
            }

            // Boundaries checks. If we hit a wall, "fake" the prev_pos such that
            // it is reflected beyond the boundary. This is done through some tricky
            // math, i.e. starting with the base velocity formulas where x is the
            // distance in one direction:
            //     velocity = (x - xₙ₋₁) / dt
            // now we want to find the fake previous pos which would be from the
            // fake "reflected" velocity, i.e. we wanna find x'ₙ₋₁
            //     reflected_velocity = (x' - x'ₙ₋₁) / dt
            // the new position would be the boundary since this is a bounce:
            //     reflected_velocity = (boundary - x'ₙ₋₁) / dt
            // so solving for prev_pos:
            //     x'ₙ₋₁ = boundary - reflected_velocity * dt
            // since reflected_velocity is really just -velocity, rewrite as:
            //     x'ₙ₋₁ = boundary + velocity * dt
            // plug the other side of the original velocity equation:
            //     x'ₙ₋₁ = boundary + ((x - xₙ₋₁) / dt) * dt
            // simplify...
            //     x'ₙ₋₁ = boundary + (x - xₙ₋₁)
            // finally apply some damping to the velocity:
            //     x'ₙ₋₁ = boundary + (x - xₙ₋₁) * bounciness
            let screen_width = screen_width();
            let screen_height = screen_height();
            if particle.pos.x < 0.0 {
                particle.pos.x = 0.0;
                particle.prev_pos.x =
                    (particle.pos.x - particle.prev_pos.x) * BLOB_BOUNCINESS;
            } else if particle.pos.x > screen_width {
                particle.pos.x = screen_width;
                particle.prev_pos.x = screen_width
                    + (particle.pos.x - particle.prev_pos.x) * BLOB_BOUNCINESS;
            }

            if particle.pos.y < 0.0 {
                particle.pos.y = 0.0;
                particle.prev_pos.y =
                    (particle.pos.y - particle.prev_pos.y) * BLOB_BOUNCINESS;
            } else if particle.pos.y > screen_height {
                particle.pos.y = screen_height;
                particle.prev_pos.y = screen_height
                    + (particle.pos.y - particle.prev_pos.y) * BLOB_BOUNCINESS;
            }
        }

        // Check all particle pairs for collisions and "bump" them apart
        for i in 0..self.particles.len() {
            for j in (i + 1)..self.particles.len() {
                let distance =
                    (self.particles[i].pos - self.particles[j].pos).length();
                let min_distance = BLOB_PARTICLE_RADIUS * 2.0;

                if distance < min_distance && distance > 0.0 {
                    let overlap = min_distance - distance;
                    let direction = (self.particles[i].pos
                        - self.particles[j].pos)
                        / distance;

                    let separation = direction * (overlap * 0.5);

                    self.particles[i].pos += separation;
                    self.particles[j].pos -= separation;
                }
            }
        }
    }

    pub fn get_center_pos(&self) -> Vec2 {
        return self.particles[self.center_particle_index].pos;
    }

    pub fn move_blob(&mut self, force_vec: Vec2) {
        // Since we're using Verlet integration to apply force we
        // "fake" it by moving the prev_pos to be further away
        for particle in &mut self.particles {
            particle.prev_pos -= force_vec;
        }
    }

    pub fn draw(&self) {
        // Chaikin subdivision
        let mut points: Vec<Vec2> = self
            .outline_particles_indices
            .iter()
            .map(|&i| self.particles[i].pos)
            .collect();

        for _ in 0..BLOB_SMOOTHING_PASSES {
            let mut new_points: Vec<Vec2> =
                Vec::with_capacity(points.len() * 2);
            for i in 0..points.len() {
                let a = points[i];
                let b = points[(i + 1) % points.len()];
                new_points.push(a * 0.75 + b * 0.25);
                new_points.push(a * 0.25 + b * 0.75);
            }
            points = new_points;
        }

        // Draw the smooth outline
        for i in 0..points.len() {
            let curr = points[i];
            let next = points[(i + 1) % points.len()];

            draw_circle(curr.x, curr.y, BLOB_OUTLINE_THICKNESS / 2.0, BLACK);

            draw_line(
                curr.x,
                curr.y,
                next.x,
                next.y,
                BLOB_OUTLINE_THICKNESS,
                BLACK,
            );
        }

        if DEBUG {
            for _ in 0..self.outline_particles_indices.len() {
                // Draw all springs as thin lines
                for spring in &self.springs {
                    let particle_a_pos = self.particles[spring.particle_a].pos;
                    let particle_b_pos = self.particles[spring.particle_b].pos;

                    draw_line(
                        particle_a_pos.x,
                        particle_a_pos.y,
                        particle_b_pos.x,
                        particle_b_pos.y,
                        1.0, // Thin line
                        GREEN,
                    );
                }

                // Draw all particles as small circles
                for (i, particle) in self.particles.iter().enumerate() {
                    let color = if self.outline_particles_indices.contains(&i) {
                        RED // Outline particles in red
                    } else {
                        BLUE // Internal particles in blue
                    };

                    draw_circle(
                        particle.pos.x,
                        particle.pos.y,
                        5.0, // Small radius for debug
                        color,
                    );
                }
            }
        }
    }
}
