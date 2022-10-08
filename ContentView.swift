//
//  ContentView.swift
//  iOSChirp
//
//  Created by Bolin Wu on 10/1/22.
//

import SwiftUI
import Charts


struct ContentView: View {

    
    var body: some View {
        
        Chart(data, id: \.x_val) {
            LineMark(
                x: .value("Time", $0.x_val),
                y: .value("Waveform Value", $0.y_val)
            )
        }
        .chartXAxis{
        }/*
        .chartXAxisLabel(position: .bottom, alignment: .center) {
            Text("Time (s)")
        }*/
        .chartYAxis{
            
        }
    }
    /*
    var body: some View {
        VStack {
            Image(systemName: "globe")
                .imageScale(.large)
                .foregroundColor(.accentColor)
            Text("Hello, world!")
        }
        .padding()
    }*/
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}


