use std::fs::File;

use std::time::{Instant};

pub mod zoom {
    #[derive(Copy)]
    #[derive(Clone)]
    pub struct Sector(pub f64, pub f64, pub f64, pub f64);
    


    impl Sector {
        pub fn new(a: Point, b: Point, width: usize, height: usize) -> Sector {
            let f = Point(a.0.min(b.0), a.1.min(b.1));
            let s = Point(a.0.max(b.0), a.1.max(b.1));
            Sector(f.0, s.0, f.1, s.1).fit_frame(width as f64, height as f64)
        }
        pub fn fit_frame(&mut self, w: f64, h: f64) -> Sector { // returns a sector that fits into the ratio w:h
            if h/w < (self.3-self.2)/(self.1-self.0) {
                let hw = (self.3-self.2)/h*w/2f64; // half of the width scaled down
                let hp = super::tools::mix(self.0, self.1, 0.5f64);
                Sector(hp - hw, hp + hw, self.2, self.3)
            } else {
                let hh = (self.1-self.0)/w*h/2f64; // half of the height scaled down
                let hp = super::tools::mix(self.2, self.3, 0.5f64);
                Sector(self.0, self.1, hp + hh, hp - hh)
            }
        }
        pub fn mix(&self, other: Sector, t: f64) -> Sector {
            Sector(super::tools::mix(self.0, other.0, t), 
                super::tools::mix(self.1, other.1, t), 
                super::tools::mix(self.2, other.2, t), 
                super::tools::mix(self.3, other.3, t)
            )
        }
        pub fn w(&self) -> f64 {
            (self.0 - self.1).abs()
        }
        pub fn h(&self) -> f64 {
            (self.2 - self.3).abs()
        }
    }

    #[derive(Copy)]
    #[derive(Clone)]
    pub struct Point(pub f64, pub f64);
    
    pub struct SectorZoom {
        pub from: Sector,
        pub to: Sector,
    }

    impl SectorZoom {
        pub fn new(from: Sector, to: Sector) -> SectorZoom {
            SectorZoom{from, to}
        }

        pub fn point(from: Sector, point: Point, zoom_factor: f64) -> SectorZoom {
            let z = 1f64 / zoom_factor;
            let w = z * from.w()*0.5f64;
            let h = z * from.h()*0.5f64;
            let to = Sector(point.0 - w, point.0 + w, point.1 - h, point.1 + h);

            SectorZoom{from, to}
        }

        pub fn zoom(&self, t: f64) -> Sector {
            
            let h = 1f64 - (self.to.1 - self.to.0) / (self.from.1 - self.from.0);
            let t = (1f64 - (t * (1f64 - h).log2()).exp2()) / h; // explanation below
            /*
            We want a function such that if we with x linearly from 0 to 1, we get values between 0 and 1, 
            meaning they are in the 1x1 square that is located at 0,0;

                   --- ---  1-|-------------------------%%%|                               
                    |  b|     |                     %%%%   |                               
                    |   |     |                  %%%       |                 ############# 
                    |  ---  h-|---------------%%%----------|---------######                
                    |         |            %%%             |  #####     |                  
                   a|         |          %%            #####            |                  
                    |         |        %%         ####     |            |                  
                    |         |      %%      ###           |            |                  
                    |         |    %%    ###               |            |                  
                    |         |  %%  ##                    |            |                  
                    |         | %##                        |            |                  
                   ___      0_|%___________________________|____________|__________________
                                                            1         f_inv(h)

            a: length of 1st sector normalized (is just 1, because of division by itself)
            b: length of 2nd sector normalized by length of 1st sector

            #: f(x)     = 1 - 2^(-x)
            f_inv(x) = -log2(1 - x)

            %: f_scaled(x) = (1 - 2^(-x * log2(1 - x))) / x    <-- this is the same as: (1f64 - (x * (1f64 - h).log2()).exp2()) / h
                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ is inside the 1x1 square of the graph
            */
            self.from.mix(self.to, t) // interpolating the corners of the sectors between the first and last sector    
        }
    }
}



#[inline(always)]
pub fn mandelbrot(x: f64, y: f64) -> f64 {
    // sections from this code are from Inigo Quilez: https://www.shadertoy.com/view/MltXz2
    let c2 = x*x + y*y;
    if 256f64*c2*c2 - 96f64*c2 + 32f64*x - 3f64 < 0f64 {
        return 1f64;
    }
    if 16f64*(c2 + 2f64*x + 1f64) - 1f64 < 0f64 {
        return 1f64;
    }
    
    let max_i = 2048;
    let mut i = 0;
    let (mut zr, mut zi) = (x, y);
    let (cr, ci, mut zr2, mut zi2, mut _t) = (zr, zi, 0f64, 0f64, 0f64);
    // unoptimized loop - somehow it is faster than the optimized one... no idea why
    // while zr*zr+zi*zi <= 65535f64 && i < max_i {
    //     // calculation z = z^2 + c
    //     // (zr + zi)^2 = zr^2 + 2*zr*zi - zi^2
    //     _t = zr;
    //     zr = zr*zr - zi*zi;
    //     zi = 2f64 * _t * zi;
    //     zr += cr;
    //     zi += ci;

    //     // boundary check abs(z^2) < 256^2-1.
    //     i += 1;
    // }
    
    
    
    // optimized loop from: https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set#Optimized_escape_time_algorithms    
    while zi2+zr2<=64f64 && i < max_i {
        zr  = zr2-zi2+cr;
        zi  = _t-zr2-zi2+ci;
        zr2 = zr*zr;
        zi2 = zi*zi;
        _t  = (zr+zi)*(zr+zi);
        i  += 1;
    }

    if i == max_i {
       return 1f64;
    }
    // stepped values (just the normalized iteration count)
    //return (i as f64 / max_i as f64).sqrt();

    // smooth values (by Inigo Quilez)
    // unoptimized:
    // smooth = i - log(log(length(z))/log(2))/log(2)
    //return (i as f64 - ((zr*zr+zi*zi).sqrt().ln()/2f64.ln()).ln()/2f64.ln())/max_i as f64;
    // optimized:
    return ((i as f64 - (zr2+zi2).log2().log2() + 4f64) / max_i as f64).sqrt();
}

#[inline(always)]
pub fn _mandelbrot_trap(x: f64, y: f64, trap: &dyn Fn(f64, f64) -> f64) -> f64 {
    // sections from this code are from Inigo Quilez: https://www.shadertoy.com/view/MltXz2
    let c2 = x*x + y*y;
    if 256f64*c2*c2 - 96f64*c2 + 32f64*x - 3f64 < 0f64 {
        return 1f64;
    }
    if 16f64*(c2 + 2f64*x + 1f64) - 1f64 < 0f64 {
        return 1f64;
    }
    
    let max_i = 2048;
    let mut i = 0;
    let (mut zr, mut zi) = (x, y);
    let (cr, ci, mut zr2, mut zi2, mut _t) = (zr, zi, 0f64, 0f64, 0f64);
    // unoptimized loop - somehow it is faster than the optimized one... no idea why
    // while zr*zr+zi*zi <= 65535f64 && i < max_i {
    //     // calculation z = z^2 + c
    //     // (zr + zi)^2 = zr^2 + 2*zr*zi - zi^2
    //     _t = zr;
    //     zr = zr*zr - zi*zi;
    //     zi = 2f64 * _t * zi;
    //     zr += cr;
    //     zi += ci;

    //     // boundary check abs(z^2) < 256^2-1.
    //     i += 1;
    // }
    
    
    
    // optimized loop from: https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set#Optimized_escape_time_algorithms    
    while zi2+zr2<=64f64 && i < max_i {
        zr  = zr2-zi2+cr;
        zi  = _t-zr2-zi2+ci;
        zr2 = zr*zr;
        zi2 = zi*zi;
        _t  = (zr+zi)*(zr+zi);
        i  += 1;
    }

    if i == max_i {
       return 1f64;
    }
    // stepped values (just the normalized iteration count)
    //return (i as f64 / max_i as f64).sqrt();

    // smooth values (by Inigo Quilez)
    // unoptimized:
    // smooth = i - log(log(length(z))/log(2))/log(2)
    //return (i as f64 - ((zr*zr+zi*zi).sqrt().ln()/2f64.ln()).ln()/2f64.ln())/max_i as f64;
    // optimized:
    return ((i as f64 - (zr2+zi2).log2().log2() + 4f64) / max_i as f64).sqrt();
}

pub fn mandelbrot_its(x: f64, y: f64, max_iter: u32) -> u32 {
    // sections from this code are from Inigo Quilez: https://www.shadertoy.com/view/MltXz2
    let c2 = x*x + y*y;
    if 256f64*c2*c2 - 96f64*c2 + 32f64*x - 3f64 < 0f64 {
        return max_iter-1;
    }
    if 16f64*(c2 + 2f64*x + 1f64) - 1f64 < 0f64 {
        return max_iter-1;
    }
    
    let mut i = 0;
    let (mut zr, mut zi) = (x, y);
    let (cr, ci, mut zr2, mut zi2, mut _t) = (zr, zi, 0f64, 0f64, 0f64);
    // unoptimized loop - somehow it is faster than the optimized one... no idea why
    // while zr*zr+zi*zi <= 65535f64 && i < max_i {
    //     // calculation z = z^2 + c
    //     // (zr + zi)^2 = zr^2 + 2*zr*zi - zi^2
    //     _t = zr;
    //     zr = zr*zr - zi*zi;
    //     zi = 2f64 * _t * zi;
    //     zr += cr;
    //     zi += ci;

    //     // boundary check abs(z^2) < 256^2-1.
    //     i += 1;
    // }
    
    
    
    // optimized loop from: https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set#Optimized_escape_time_algorithms    
    while zi2+zr2<=64f64 && i < max_iter {
        zr  = zr2-zi2+cr;
        zi  = _t-zr2-zi2+ci;
        zr2 = zr*zr;
        zi2 = zi*zi;
        _t  = (zr+zi)*(zr+zi);
        i  += 1;
    }
    i-1
}

pub fn render_histogram(sector: zoom::Sector, max_iter: u32, fractal: &dyn Fn(f64, f64, u32) -> u32, detail: usize) {
    let width = (detail as f64 * sector.w()) as usize;
    let height = (detail as f64 * sector.h()) as usize;
    
    let mut histogram = vec![0u32; max_iter as usize];

    for x in 0..(width) as usize {
        for y in 0..(height) as usize {
            histogram[fractal(x as f64 / detail as f64 + sector.0, y as f64 / detail as f64 + sector.2, max_iter) as usize] += 1;
        }
    }
    let mut max = 0u32;
    for i in histogram.iter() {
        max = *i.max(&max)
    }
    max = (max as f64/2048f64) as u32;
    let mut image = vec![0u8; (max_iter*max*4) as usize];
    for x in 0..max_iter as usize {
        let h = histogram[x];
        for y in 0..max as usize {
            let col = if max as usize - y <= (h as f64 / 2048f64) as usize {1} else {0};
            for i in 0..4 {
                image[i+(x+y*max_iter as usize)*4] = col*255;
            }
        }
    }
    save_frame(&image, max_iter as usize, max as usize, 0);
}

pub fn render_frame(
        col: &mut Vec<u8>, 
        width: usize, height: usize, 
        fractal: &dyn Fn(f64, f64) -> f64, 
        color: &dyn Fn(f64) -> [f64; 4], 
        sector: zoom::Sector,
    ) {
    if width*height*4 != col.len() {
        panic!("color buffer needs to be 4 times as big as the width and height multiplied (because of RGBA mode)")
    }
    
    let mut values = vec![0f64; width*height];
    let (w, h) = (width as f64, height as f64);
    for (i, val) in values.iter_mut().enumerate() {
        let x = (i as f64 / w).fract();
        let y = (i as f64 / w).floor() / h;
        *val = fractal(tools::mix(sector.0, sector.1, x), tools::mix(sector.2, sector.3, y));
    }
    for (i, col) in col.chunks_exact_mut(4).enumerate() {
        let v = color(*values.get(i).unwrap());
        col[0] = (v[0] * 255f64) as u8;
        col[1] = (v[1] * 255f64) as u8;
        col[2] = (v[2] * 255f64) as u8;
        col[3] = (v[3] * 255f64) as u8;
    }
}

pub fn render_sequence(zoom: zoom::SectorZoom, width: usize, height: usize, frames: u32) {
    
    let mut color_data = vec![0u8; width*height*4];

    if frames == 0 {
        render_frame(&mut color_data, width, height, &mandelbrot, &palette, zoom.to);
        save_frame(&color_data, width, height, 0);
        return
    }

    for frame_count in 0..frames+1 {
        println!("frame: {:3}/{:3}", frame_count, frames);

        let mut timer = Instant::now();
        render_frame(&mut color_data, width, height, &mandelbrot, &palette, zoom.zoom(frame_count as f64 / frames as f64));
        
        println!("te mandelbrot: {:06} ms", timer.elapsed().as_millis());

        
        timer = Instant::now();

        save_frame(&color_data, width, height, frame_count);
        
        println!("te saving:     {:06} ms\n", timer.elapsed().as_millis());
    }
}

pub fn save_frame(buf: &Vec<u8>, width: usize, height: usize, index: u32) {
    let str_path = format!("output/pic{:06}.png", index);
        
    let path = std::path::Path::new(&str_path);
    let file = File::create(path).unwrap();
    let ref mut w = std::io::BufWriter::new(file);

    let mut encoder = png::Encoder::new(w, width as u32, height as u32);
    encoder.set_color(png::ColorType::RGBA);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header().unwrap();

    writer.write_image_data(&buf[..]).unwrap();
}

pub fn palette(x: f64) -> [f64; 4] {
    //from here: https://www.shadertoy.com/view/tsKcWm
    let t = x.fract()*2f64;
    let x = 1f64 - (1f64 - t).abs();
    let sx = x.sqrt();
    let xx = x * x;
    if t < 1f64 {
        return [sx, x, xx, 1f64];
    } else {
        return [xx, x, sx, 1f64];
    }
}

pub mod tools {
    pub fn mix(a: f64, b: f64, t: f64) -> f64 {
        b * t + a * (1f64 - t)
    }
    
    pub fn smooth(x: f64) -> f64 {
        x*x*(3f64-2f64*x)
    }
}
