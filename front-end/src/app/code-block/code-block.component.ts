import { Component, Input } from '@angular/core';

import { Block } from '../block.interface';
import { BlockService } from '../block.service';

@Component({
  selector: 'app-code-block',
  templateUrl: './code-block.component.html',
  styleUrls: ['./code-block.component.css'],
})

export class CodeBlockComponent {
  @Input() block!: Block;
  executingBlocks: boolean = false;

  constructor(private blockService: BlockService) { 
    this.blockService.executingBlocks.subscribe(
      (res) => { this.executingBlocks = res; },
    );
  }

  removeBlock(): void {
    this.blockService.removeBlock(this.block.blockId);
  }
}
