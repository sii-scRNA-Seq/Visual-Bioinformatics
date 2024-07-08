import { Component } from '@angular/core';

import { Block, BlockId } from '../block.interface';
import { BlockService } from '../block.service';
import { OutputService } from '../output.service';

@Component({
  selector: 'app-block-library',
  templateUrl: './block-library.component.html',
  styleUrls: ['./block-library.component.css'],
})

export class BlockLibraryComponent {
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

  addBlock(id: BlockId): void {
    this.blockService.addBlock(id);
  }
}
