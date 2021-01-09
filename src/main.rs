mod fractals;
use fractals::zoom::*;


fn main() {
    // to convert the frames to a video:
    // ffmpeg -framerate 24 -i output/pic%06d.png mandelbrot.mp4
    let width: usize = 1920;
    let height: usize = 1080;
    // let width: usize = 192*2;
    // let height: usize = 108*2;
    let frames = 0;
    
    let fs = Sector::new(Point(-2f64, -2f64), Point(2f64, 2f64), width, height);                       // sector of first frame
    //let ls = Sector::new(Point(0.2550008f64, -0.0006f64), Point(0.2550009f64, -0.000599f64), width, height); // sector of last frame
    //let zoom = SectorZoom::new(fs, ls);
    let zoom = SectorZoom::point(fs, Point(-1.63251010322,-0.0), 100000000000f64);
    fractals::render_sequence(zoom, width, height, frames);
    // fractals::render_histogram(fs, 1024, &fractals::mandelbrot_its, 1024);
}
