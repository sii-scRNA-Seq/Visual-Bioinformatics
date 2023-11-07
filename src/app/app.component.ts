import { Component, OnInit } from '@angular/core';
import { OutputService } from './output.service';
import { Output } from './output'

@Component({
  selector: 'app-root',
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.css']
})

export class AppComponent implements OnInit {

  text: String = "...";

  constructor(private outputService:OutputService) { }

  ngOnInit() {
	  this.getOutputs();
  }

  getOutputs() {
	  this.outputService.getOutputs().subscribe((data: Output) => this.text = data.title);
  }
  
}