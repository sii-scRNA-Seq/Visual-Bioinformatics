import { Component } from '@angular/core';

import { BlockId } from '../block.interface';
import { BlockService } from '../block.service';

@Component({
  selector: 'app-block-library',
  templateUrl: './block-library.component.html',
  styleUrls: ['./block-library.component.css'],
})

export class BlockLibraryComponent {
  constructor(private blockService: BlockService) { }

  addBlock(id: BlockId): void {
    this.blockService.addBlock(id);
  }
}
