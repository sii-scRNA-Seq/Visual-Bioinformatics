import { Component } from '@angular/core';

import { BlockIdToTitleMap } from '../block.interface';
import { Output } from './../output';
import { OutputService } from '../output.service';

@Component({
  selector: 'app-output-display',
  templateUrl: './output-display.component.html',
  styleUrls: ['./output-display.component.css'],
})

export class OutputDisplayComponent {

  BlockIdToTitleMap = BlockIdToTitleMap;

  outputList: Output[] = [];

  constructor(private outputService: OutputService) {
    this.outputService.outputs.subscribe(
      (res) => { this.outputList = res; }
    );
  }
}
