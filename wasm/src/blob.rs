// Blob constants
const BLOB_SPRING_STIFFNESS: f32 = 40.0;
const BLOB_SHAPE_STIFFNESS: f32 = 96.0;
const BLOB_BOUNCINESS: f32 = 0.1; // Sane values are 0.0 - 1.0
const BLOB_RADIUS: f32 = 100.0;
const BLOB_PARTICLE_RADIUS: f32 = 16.0;
const BLOB_MASS: f32 = 32.0;
const BLOB_OUTLINE_THICKNESS: f32 = 24.0;
const GRAVITY: f32 = 1600.0;
const BLOB_SPRING_DAMPING: f32 = BLOB_SPRING_STIFFNESS / 1.5;
const VELOCITY_DAMPING: f32 = 0.99;
const BLOB_MAX_SPEED: f32 = 160.0;
const EPSILON: f32 = 0.00000001;
use macroquad::{
    color::{BLACK, BLUE, GREEN, PINK, RED},
    math::Vec2,
    rand,
    shapes::{draw_circle, draw_line},
    window::{screen_height, screen_width},
};
use std::f32::consts::{PI, SQRT_2};

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
        // Poisson disk sampling to pack particles
        const SAMPLES: usize = 100;
        let cell_size = BLOB_PARTICLE_RADIUS * SQRT_2.recip();
        let grid_width = (BLOB_RADIUS * 2.0 / cell_size).ceil() as usize;
        let grid_height = grid_width;

        // For checking if there are samples too close
        let mut grid: Vec<Vec<Option<Vec2>>> =
            vec![vec![None; grid_width]; grid_height];

        // Will hold "active points"
        let mut active_list: Vec<Vec2> = Vec::new();

        // Pick the first sample in the center and add to the queue
        let first_sample = Vec2::new(0.0, 0.0);
        let grid_x = ((first_sample.x + BLOB_RADIUS) / cell_size) as usize;
        let grid_y = ((first_sample.y + BLOB_RADIUS) / cell_size) as usize;
        grid[grid_x][grid_y] = Some(first_sample);
        active_list.push(first_sample);

        while !active_list.is_empty() {
            let i =
                (rand::gen_range(0.0, 1.0) * active_list.len() as f32) as usize;
            let parent = active_list[i];

            // Try to generate k candidates around this parent
            let mut found = false;
            for j in 0..SAMPLES {
                // Generate candidates at random angles, just far enough away
                let angle = 2.0
                    * PI
                    * (rand::gen_range(0.0, 1.0) + j as f32 / SAMPLES as f32);
                let radius = BLOB_PARTICLE_RADIUS * 2.0 + EPSILON;
                let x = parent.x + radius * angle.cos();
                let y = parent.y + radius * angle.sin();
                let candidate = Vec2 { x: x, y: y };
                let distance_from_center = candidate.length();
                if distance_from_center <= BLOB_RADIUS - radius {
                    // Check if candidate is far enough from existing samples
                    let candidate_grid_x =
                        ((candidate.x + BLOB_RADIUS) / cell_size) as usize;
                    let candidate_grid_y =
                        ((candidate.y + BLOB_RADIUS) / cell_size) as usize;

                    let mut is_far_enough = true;

                    'outer: for step_x in -2..=2 {
                        let check_x = candidate_grid_x as i32 + step_x;
                        if check_x < 0 || check_x >= grid_width as i32 {
                            // Out of X bounds
                            continue;
                        }

                        for step_y in -2..=2 {
                            let check_y = candidate_grid_y as i32 + step_y;
                            if check_y < 0 || check_y >= grid_height as i32 {
                                // Out of Y bounds
                                continue;
                            }

                            if let Some(existing_sample) =
                                grid[check_x as usize][check_y as usize]
                            {
                                let distance =
                                    (candidate - existing_sample).length();
                                if distance < BLOB_PARTICLE_RADIUS * 2.0 {
                                    is_far_enough = false;
                                    break 'outer;
                                }
                            }
                        }
                    }

                    if is_far_enough {
                        // Found a valid candidate, add to grid and active_list
                        grid[candidate_grid_x][candidate_grid_y] =
                            Some(candidate);
                        active_list.push(candidate);
                        found = true;
                        break; // Break from SAMPLES loop
                    }
                }
            }

            if !found {
                active_list.swap_remove(i);
            }
        }

        // Create particles pased on the Poisson disc sampling
        let inner_particle_positions: Vec<Vec2> =
            grid.iter().flatten().filter_map(|&sample| sample).collect();

        let num_inner_particles = inner_particle_positions.len();

        // Find the "center" particle (the one with min distance from origin)
        let mut center_particle_index = 0;
        let mut min_distance = (inner_particle_positions[0] - origin).length();

        for (i, &pos) in inner_particle_positions.iter().enumerate().skip(1) {
            let distance = pos.length();
            if distance < min_distance {
                min_distance = distance;
                center_particle_index = i;
            }
        }

        // Add outline particles
        let circumference = 2.0 * PI * BLOB_RADIUS;
        let num_outline_particles =
            (circumference / (BLOB_PARTICLE_RADIUS * 2.0)).round() as usize;

        let outline_particle_positions: Vec<Vec2> = (0..num_outline_particles)
            .map(|i| {
                let angle = 2.0 * PI * i as f32 / num_outline_particles as f32;
                Vec2::new(BLOB_RADIUS * angle.cos(), BLOB_RADIUS * angle.sin())
            })
            .collect();

        let outline_particle_indices: Vec<usize> = (num_inner_particles
            ..num_inner_particles + num_outline_particles)
            .collect();

        let all_particle_positions: Vec<Vec2> =
            [inner_particle_positions, outline_particle_positions].concat();

        let center_of_mass = all_particle_positions.iter().sum::<Vec2>()
            / all_particle_positions.len() as f32;

        // Create all the particles, translating to the origin
        let particles: Vec<Particle> = all_particle_positions
            .iter()
            .map(|&pos| Particle {
                pos: origin + pos,
                prev_pos: origin + pos,
                offset_from_center: pos - center_of_mass,
            })
            .collect();

        // Create springs connecting particles to the closest 6 inner particles
        let mut springs: Vec<Spring> = Vec::new();
        for i in 0..particles.len() {
            let mut candidates: Vec<Spring> = (0..particles.len())
                .filter(|&j| j != i && !outline_particle_indices.contains(&j))
                .map(|j| Spring {
                    particle_a: i,
                    particle_b: j,
                    rest_length: (particles[i].pos - particles[j].pos).length(),
                })
                .collect();

            candidates.sort_by(|spring_a, spring_b| {
                spring_a
                    .rest_length
                    .partial_cmp(&spring_b.rest_length)
                    .unwrap()
            });

            let num_connections = if outline_particle_indices.contains(&i) {
                6
            } else {
                8
            };
            springs.extend(candidates.into_iter().take(num_connections));
        }

        // Connect outline particles to their neighbors
        for i in 0..outline_particle_indices.len() {
            let current_idx = outline_particle_indices[i];
            let next_idx = outline_particle_indices
                [(i + 1) % outline_particle_indices.len()];

            springs.push(Spring {
                particle_a: current_idx,
                particle_b: next_idx,
                rest_length: (particles[current_idx].pos
                    - particles[next_idx].pos)
                    .length(),
            });
        }

        Blob {
            particles: particles,
            springs: springs,
            outline_particles_indices: outline_particle_indices,
            center_particle_index: center_particle_index,
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
        for i in 0..self.outline_particles_indices.len() {
            let prev = (i + self.outline_particles_indices.len() - 1)
                % self.outline_particles_indices.len();
            let curr = i;
            let next = (i + 1) % self.outline_particles_indices.len();

            let pos_prev =
                self.particles[self.outline_particles_indices[prev]].pos;
            let pos_curr =
                self.particles[self.outline_particles_indices[curr]].pos;
            let pos_next =
                self.particles[self.outline_particles_indices[next]].pos;

            // Smooth current position by averaging with neighbors
            let smooth_pos = (pos_prev + pos_curr * 2.0 + pos_next) * 0.25;

            let next_smooth = (pos_curr
                + pos_next * 2.0
                + self.particles[self.outline_particles_indices
                    [(next + 1) % self.outline_particles_indices.len()]]
                .pos)
                * 0.25;

            draw_circle(
                smooth_pos.x,
                smooth_pos.y,
                BLOB_OUTLINE_THICKNESS / 2.0,
                BLACK,
            );

            draw_line(
                smooth_pos.x,
                smooth_pos.y,
                next_smooth.x,
                next_smooth.y,
                BLOB_OUTLINE_THICKNESS,
                BLACK,
            );

            if DEBUG {
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
                    } else if i == self.center_particle_index {
                        PINK
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
