import { Component, OnInit } from '@angular/core';
import { OutputService } from './../output.service';
import { Output } from './../output';

@Component({
  selector: 'app-output-display',
  templateUrl: './output-display.component.html',
  styleUrls: ['./output-display.component.css'],
})

export class OutputDisplayComponent implements OnInit {
  text = '...';

  constructor(private outputService: OutputService) { }

  ngOnInit() {
    this.getOutputs();
  }

  getOutputs() {
    this.outputService.getOutputs().subscribe((data: Output) => this.text = data.title);
  }
}
