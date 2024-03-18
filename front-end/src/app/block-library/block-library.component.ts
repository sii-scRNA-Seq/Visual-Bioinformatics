import { Component } from '@angular/core';

import { BlockId } from '../block.interface';
import { BlockService } from '../block.service';

@Component({
  selector: 'app-block-library',
  templateUrl: './block-library.component.html',
  styleUrls: ['./block-library.component.css'],
})

export class BlockLibraryComponent {
  executingBlocks: boolean = false;

  constructor(private blockService: BlockService) {
    this.blockService.executingBlocks.subscribe(
      (res) => { this.executingBlocks = res; },
    );
  }

  addBlock(id: BlockId): void {
    this.blockService.addBlock(id);
  }
}
