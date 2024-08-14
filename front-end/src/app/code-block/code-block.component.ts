import { Component, Input } from '@angular/core';

import { Block } from '../block.interface';
import { BlockService } from '../block.service';
import { CurrentDatasetService } from '../current-dataset.service';
import { OutputService } from '../output.service';

@Component({
  selector: 'app-code-block',
  templateUrl: './code-block.component.html',
  styleUrls: ['./code-block.component.css'],
})

export class CodeBlockComponent {
  @Input() block!: Block;
  executingBlocks: boolean = false;

  constructor(private blockService: BlockService, private currentDatasetService: CurrentDatasetService, private outputService: OutputService) { 
    this.outputService.executingBlocks.subscribe(
      (res) => { this.executingBlocks = res; },
    );
  }

  onMatSelectValueChanges(): void {
    this.currentDatasetService.onMatSelectValueChanges(this.block);
  }

  removeBlock(): void {
    this.blockService.removeBlock(this.block.blockUUID);
  }
}
