import { Component } from '@angular/core';
import { OutputService } from './../output.service';
import { Output } from './../output';

@Component({
  selector: 'app-output-display',
  templateUrl: './output-display.component.html',
  styleUrls: ['./output-display.component.css'],
})

export class OutputDisplayComponent {
  outputs: Output[] = [];

  constructor(private outputService: OutputService) {
    this.outputService.outputs.subscribe(
      (res) => { this.outputs = res; },
    );
  }
}
