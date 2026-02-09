use macroquad::prelude::*;
mod blob;

#[macroquad::main("blob")]
async fn main() {
    let initial_screen_width = screen_width();
    let initial_screen_height = screen_height();

    let mut blob = blob::Blob::new(Vec2 {
        x: initial_screen_width / 2.0,
        y: initial_screen_height / 2.0,
    });

    let mut time: f32 = 0.0;

    loop {
        let dt = get_frame_time().min(1.0 / 120.0);
        time += dt;

        clear_background(WHITE);
        draw_text(&format!("FPS: {}", get_fps()), 0., 16., 32., BLACK);

        if is_mouse_button_down(MouseButton::Left) {
            let mouse_pos = Vec2::from(mouse_position());
            let blob_pos = blob.get_center_pos();
            let direction = (mouse_pos - blob_pos) * dt * 0.20; // Scale it down
            blob.move_blob(direction);

            draw_poly_lines(
                mouse_pos.x,
                mouse_pos.y,
                3,
                18.0,
                time * 200.0,
                20.0,
                RED,
            );
        }

        blob.update(get_frame_time().min(1.0 / 120.0));
        blob.draw();
        next_frame().await
    }
}
