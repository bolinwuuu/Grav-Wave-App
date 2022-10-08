//
//  run_chirp.swift
//  iOSChirp
//
//  Created by Bolin Wu on 10/1/22.
//

import Foundation
import UIKit
import Accelerate
import Darwin
import SwiftUI
import Charts


/*
extension Array {
    func scan<T>(initial: T, _ f: (T, Element) -> T) -> [T] {
        return self.reduce([initial], { (listSoFar: [T], next: Element) -> [T] in
            // because we seeded it with a non-empty
            // list, it's easy to prove inductively
            // that this unwrapping can't fail
            let lastElement = listSoFar.last!
            return listSoFar + [f(lastElement, next)]
        })
    }
}
extension Sequence {
    func scan<Result>(
        _ initial: Result,
        _ nextPartialResult: (Result, Element) -> Result
    ) -> [Result] {
        var result = [initial]
        for x in self {
            result.append(nextPartialResult(result.last!, x))
        }
        return result
    }
    
    
}


func cumsum(inarr: [Float]) -> [Float] {
    var outarr = [Float](repeating: 0, count: inarr.count);
    outarr[0] = inarr[0];
    for index in 1...inarr.count - 1 {
        outarr[index] = outarr[index - 1] + inarr[index];
    }
    return outarr;
}*/

struct MonthlyHoursOfSunshine {
    var date: Date
    var hoursOfSunshine: Double

    init(month: Int, hoursOfSunshine: Double) {
        let calendar = Calendar.autoupdatingCurrent
        self.date = calendar.date(from: DateComponents(year: 2020, month: month))!
        self.hoursOfSunshine = hoursOfSunshine
    }
}


struct Coords {
    var x_val: Double
    var y_val: Double
    
    init(x_in: Double, y_in: Double) {
        self.x_val = x_in;
        self.y_val = y_in;
    }
}


class Run_Chirp {
    
    var m1: Double
    var m2: Double
    
    var sampleN: Double
    var t: [Double]
    var freq: [Double]
    var h: [Double]

    init(mass1: Double, mass2: Double) {
        m1 = mass1;
        m2 = mass2;
        
        // Implied chirp mass (governs frequency and amplitude evolution)
        // (PPNP text right after Eqn 74)
          
        var mchirp:Double = pow((m1*m2),(3/5))/pow((m1+m2),(1/5));

        // Physical constants

        let g:Double = 6.67e-11;
        let c:Double = 2.998e8;
        let pc:Double = 3.086e16;
        let msun:Double = 2.0e30;

        // Compute Schwarzchild radii of stars

        var r1 = 2 * g * m1 * msun / pow(c, 2);
        var r2 = 2 * g * m2 * msun / pow(c, 2);

        // Frequency coefficient
        // (Based on PPNP Eqn 73)

        var fcoeff:Double = (1/(8*Double.pi)) * pow(pow(5, 3), 1/8) * pow(pow(c, 3) / (g*mchirp*msun), 5/8);

        // Amplitude coefficient (assume source at 15 Mpc)
        // (Based on PPNP Eqn 74)

        var rMpc:Double = 15;
        var r = rMpc * 1e6 * pc;
        var hcoeff = (1/r) * pow(5*pow(g*mchirp*msun/pow(c, 2), 5)/c, 1/4);

        // Amplitude rescaling parameter

        let hscale:Double = 1e21;

        // frequency (Hz) when signal enters detector band

        let fbandlo:Double = 30;

        // Compute time remaining to coalescence from entering band
        // (Based on PPNP Eqn 73)

        var tau = pow(fcoeff/fbandlo, 8/3);

        // Debugging summary

        print("Starting chirp simulation with M1, M2, Mchirp = " + String(m1) + " " + String(m2) + " " + String(mchirp) + " " + "(Msun)");
        print("--> Schwarzchild radii = " + String(r1) + " " + String(r2) + "m");
        print("Distance to source r = " + String(rMpc) + " Mpc");
        print("Detection band low frequency = " + String(fbandlo) + "Hz\n--> Time to coalescence = " + String(tau) + " s\n");

        // Sampling rate (Hz) - fixed at 48 kHz for mp4 output

        let fsamp:Double = 48000;
        var dt = 1/fsamp;

        // Length of time to simulate (round up to nearest tenth of an integer second and add a tenth)

        var upperT = ceil(10*tau)/10 + 0.1;

        // Create time sample container

        sampleN = floor(fsamp*upperT);
        //var t = Array(0...Int(upperN)-1);
        t = Array(stride(from: 0, through: sampleN-1, by: 1));
        t = vDSP.multiply(dt, t);

        // Determine frequency (and then time) when Schwarzchild radii touch
        // (Use Kepler's 3rd law)
        // (Double orbital frequency to get GW frequency)

        var ftouch = 2 * (1/(2*Double.pi)) * pow(g*(m1+m2)*msun/pow(r1+r2,3), 1/2);
        var tautouch = pow(fcoeff/ftouch, 8/3);
        print("GW frequency when Schwarzchild radii touch: " + String(ftouch) + " Hz\n--> Occurs " + String(tautouch) + " seconds before point-mass coalescence\n");
        
        // Create frequency value vs time (up to last time sample before point-mass coalescence)
        // (Based on PPNP Eqn 73)
        
        //var minusdt = -dt;
        

        var vzero:Double = 0;
        var iTau:Double = floor(tau / dt);

        var lastSample = floor((pow(ftouch / fcoeff, -8/3) - tau) / -dt);
        var maxFreq:Double = pow(-lastSample * dt + tau, -3/8) * fcoeff;
        
        var freq1 = Array(stride(from: 0, through: lastSample, by: 1));
        vDSP.multiply(-dt, freq1, result: &freq1);
        //var freq1 = [Float](repeating: 0, count: Int(iTau) + 1);
        var freq2 = [Double](repeating: maxFreq, count: Int(sampleN - lastSample) - 1);
        /*
        vDSP_vramp(&vzero,
                   &minusdt,
                   &freq1,
                   vDSP_Stride(1),
                   vDSP_Length(iTau + 1));*/

        vDSP.add(tau, freq1, result: &freq1);
        //freq = freq.map{pow(Float($0), -3/8)};
        
        var exp = [Double](repeating: -3/8, count: freq1.count);
        freq = vForce.pow(bases: freq1, exponents: exp);
        vDSP.multiply(fcoeff, freq, result: &freq);
        freq += freq2;

        //Create amplitude value vs time (up to last time sample before touch)
        // (Based on PPNP Eqn 74)
        
        exp = [Double](repeating: -1/4, count: freq1.count);
        var amp = vForce.pow(bases: freq1, exponents: exp);
        amp = vDSP.multiply(hcoeff * hscale, amp);
        var amp2 = [Double](repeating: 0, count: Int(sampleN - lastSample) - 1);
        amp += amp2;
        
        // Generate strain signal in time domain

        var phi = [Double](repeating: 0, count: freq.count);
        // Cumulative sum of freq
        phi[0] = freq[0];
        for index in 1...freq.count - 1 {
            phi[index] = phi[index - 1] + freq[index];
        }
        vDSP.multiply(2 * Double.pi * dt, phi, result: &phi);
        
        h = vDSP.multiply(amp, vForce.sin(phi));
    }
    
    func run_mode(runMode: Int) -> [Coords] {
        if (runMode == 2) {
            return run_mode_2();
        }
        else if (runMode == 3) {
            return run_mode_3();
        }
        else {
            var ret = [Coords(x_in: 0, y_in: 0)];
            return ret;
        }
    }
    
    func run_mode_2() -> [Coords] {
        var cd = [Coords](repeating: Coords(x_in: 0, y_in: 0), count: self.t.count);
        
        var idx = 0;
        while (idx < cd.count) {
            cd[idx].x_val = self.t[idx];
            cd[idx].y_val = self.h[idx];
            idx += 1;
        }
        
        return cd;
    }
    
    func run_mode_3() -> [Coords] {
        var cd = [Coords](repeating: Coords(x_in: 0, y_in: 0), count: self.t.count);
        
        var idx = 0;
        while (idx < cd.count) {
            cd[idx].x_val = self.t[idx];
            cd[idx].y_val = self.freq[idx];
            idx += 1;
        }
        
        return cd;
    }
    
    func run_mode_4() {
        vDSP.convert(amplitude: freq, toDecibels: &freq, zeroReference: Double(freq.count));
    }
    
    
};

var test = Run_Chirp(mass1: 5, mass2: 20);

var data = test.run_mode(runMode: 2);

/*
var body: some View {
    Chart(data, id: \.x_val) {
        LineMark(
            x: .value("Time", $0.x_val),
            y: .value("H", $0.y_val)
        )
    }
}*/






/*
func run_chirp(mass1:Double, mass2:Double, runMode:Int) {
    
    print("Entered run_chirp with\nMass1 = " + String(mass1) + "\nMass2 = " + String(mass2) + "\nRunMode = " + String(runMode) + "\n");

    // var option:Int = 5;
    var m1:Double = mass1;
    var m2:Double = mass2;
    /*
    if (option == 1) {
        
    // Neutron star masses:
        
        m1 = 1.4;
        m2 = 1.4;
    } else if (option == 2) {
        
    // Light black holes:
        
        m1 = 5;
        m2 = 5;
    } else if (option==3) {

    // Typical black holes:

        m1 = 10;
        m2 = 10;
    } else if (option == 4) {

    // Heavy black holes:

        m1 = 30;
        m2 = 30;
    } else if (option == 5) {

    // Customized:

        m1 = 5;
        m2 = 20;
    }*/

    // Implied chirp mass (governs frequency and amplitude evolution)
    // (PPNP text right after Eqn 74)
      
    var mchirp:Double = pow((m1*m2),(3/5))/pow((m1+m2),(1/5));

    // Physical constants

    let g:Double = 6.67e-11;
    let c:Double = 2.998e8;
    let pc:Double = 3.086e16;
    let msun:Double = 2.0e30;

    // Compute Schwarzchild radii of stars

    var r1 = 2 * g * m1 * msun / pow(c, 2);
    var r2 = 2 * g * m2 * msun / pow(c, 2);

    // Frequency coefficient
    // (Based on PPNP Eqn 73)

    var fcoeff:Double = (1/(8*Double.pi)) * pow(pow(5, 3), 1/8) * pow(pow(c, 3) / (g*mchirp*msun), 5/8);

    // Amplitude coefficient (assume source at 15 Mpc)
    // (Based on PPNP Eqn 74)

    var rMpc:Double = 15;
    var r = rMpc * 1e6 * pc;
    var hcoeff = (1/r) * pow(5*pow(g*mchirp*msun/pow(c, 2), 5)/c, 1/4);

    // Amplitude rescaling parameter

    let hscale:Double = 1e21;

    // frequency (Hz) when signal enters detector band

    let fbandlo:Double = 30;

    // Compute time remaining to coalescence from entering band
    // (Based on PPNP Eqn 73)

    var tau = pow(fcoeff/fbandlo, 8/3);

    // Debugging summary

    print("Starting chirp simulation with M1, M2, Mchirp = " + String(m1) + " " + String(m2) + " " + String(mchirp) + " " + "(Msun)");
    print("--> Schwarzchild radii = " + String(r1) + " " + String(r2) + "m");
    print("Distance to source r = " + String(rMpc) + " Mpc");
    print("Detection band low frequency = " + String(fbandlo) + "Hz\n--> Time to coalescence = " + String(tau) + " s\n");

    // Sampling rate (Hz) - fixed at 48 kHz for mp4 output

    let fsamp:Double = 48000;
    var dt = 1/fsamp;

    // Length of time to simulate (round up to nearest tenth of an integer second and add a tenth)

    var upperT = ceil(10*tau)/10 + 0.1;

    // Create time sample container

    var upperN = floor(fsamp*upperT);
    //var t = Array(0...Int(upperN)-1);
    var t = Array(stride(from: 0, through: upperN-1, by: 1));
    t = vDSP.multiply(dt, t);

    // Determine frequency (and then time) when Schwarzchild radii touch
    // (Use Kepler's 3rd law)
    // (Double orbital frequency to get GW frequency)

    var ftouch = 2 * (1/(2*Double.pi)) * pow(g*(m1+m2)*msun/pow(r1+r2,3), 1/2);
    var tautouch = pow(fcoeff/ftouch, 8/3);
    print("GW frequency when Schwarzchild radii touch: " + String(ftouch) + " Hz\n--> Occurs " + String(tautouch) + " seconds before point-mass coalescence\n");
    
    // Create frequency value vs time (up to last time sample before point-mass coalescence)
    // (Based on PPNP Eqn 73)
    
    //var minusdt = -dt;
    

    var vzero:Double = 0;
    var iTau:Double = floor(tau / dt);

    var lastSample:Double = floor((pow(ftouch / fcoeff, -8/3) - tau) / -dt);
    var maxFreq:Double = pow(-lastSample * dt + tau, -3/8) * fcoeff;
    
    var freq1 = Array(stride(from: 0, through: lastSample, by: 1));
    freq1 = vDSP.multiply(-dt, freq1);
    //var freq1 = [Float](repeating: 0, count: Int(iTau) + 1);
    var freq2 = [Double](repeating: maxFreq, count: Int(upperN - lastSample) - 1);
    /*
    vDSP_vramp(&vzero,
               &minusdt,
               &freq1,
               vDSP_Stride(1),
               vDSP_Length(iTau + 1));*/

    freq1 = vDSP.add(tau, freq1);
    //freq = freq.map{pow(Float($0), -3/8)};
    
    var exp = [Double](repeating: -3/8, count: freq1.count);
    var freq = vForce.pow(bases: freq1, exponents: exp);
    freq = vDSP.multiply(fcoeff, freq);
    freq += freq2;

    //Create amplitude value vs time (up to last time sample before touch)
    // (Based on PPNP Eqn 74)
    
    exp = [Double](repeating: -1/4, count: freq1.count);
    var amp = vForce.pow(bases: freq1, exponents: exp);
    amp = vDSP.multiply(hcoeff * hscale, amp);
    var amp2 = [Double](repeating: 0, count: Int(upperN - lastSample) - 1);
    amp += amp2;
    
    // Generate strain signal in time domain

    var phi = [Double](repeating: 0, count: freq.count);
    // Cumulative sum of freq
    phi[0] = freq[0];
    for index in 1...freq.count - 1 {
        phi[index] = phi[index - 1] + freq[index];
    }
    phi = vDSP.multiply(2 * Double.pi * dt, phi);
    
    var h = vDSP.multiply(amp, vForce.sin(phi));

    /*
    var body: some View {
        Chart() {
            LineMark(
                x: .value("Time", t),
                y: .value("Waveform Value", h)
            )
        }
    }*/

    

    
}



var data: [MonthlyHoursOfSunshine] = [
    MonthlyHoursOfSunshine(month: 1, hoursOfSunshine: 74),
    MonthlyHoursOfSunshine(month: 2, hoursOfSunshine: 99),
    MonthlyHoursOfSunshine(month: 3, hoursOfSunshine: 45),
    MonthlyHoursOfSunshine(month: 4, hoursOfSunshine: 65),
    MonthlyHoursOfSunshine(month: 5, hoursOfSunshine: 83),
    MonthlyHoursOfSunshine(month: 6, hoursOfSunshine: 44),
    MonthlyHoursOfSunshine(month: 7, hoursOfSunshine: 92),
    MonthlyHoursOfSunshine(month: 8, hoursOfSunshine: 62)
]*/



