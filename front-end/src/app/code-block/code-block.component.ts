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

  constructor(private blockService: BlockService) { }

  removeBlock(): void {
    this.blockService.removeBlock(this.block.blockId);
  }
}
