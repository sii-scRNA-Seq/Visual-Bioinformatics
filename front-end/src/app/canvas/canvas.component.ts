import { Component } from '@angular/core';

import { Block } from '../block.interface';
import { BlockService } from '../block.service';
import { OutputService } from '../output.service';

@Component({
  selector: 'app-canvas',
  templateUrl: './canvas.component.html',
  styleUrls: ['./canvas.component.css'],
})

export class CanvasComponent {
  blockList: Block[] = [];
  executingBlocks: boolean = false;

  constructor(private blockService: BlockService, private outputService: OutputService) {
    this.blockService.blocksOnCanvas.subscribe(
      (res) => { this.blockList = res; },
    );
    this.outputService.executingBlocks.subscribe(
      (res) => { this.executingBlocks = res; },
    );
  }

  executeBlocks() {
    this.blockService.executeBlocks();
  }
}
