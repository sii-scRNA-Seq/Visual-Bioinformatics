import { Component, Input } from '@angular/core';

import { Block } from '../block.interface';
import { BlockService } from '../block.service';
import { DatasetInfo } from '../dataset-info';
import { DatasetInfoService } from '../dataset-info.service';
import { OutputService } from '../output.service';

@Component({
  selector: 'app-code-block',
  templateUrl: './code-block.component.html',
  styleUrls: ['./code-block.component.css'],
})

export class CodeBlockComponent {
  @Input() block!: Block;
  datasetInfo: DatasetInfo[] = [];
  executingBlocks: boolean = false;

  constructor(private blockService: BlockService, private datasetInfoService: DatasetInfoService, private outputService: OutputService) { 
    this.datasetInfoService.datasetInfo.subscribe(
      (res) => { this.datasetInfo = res; },
    );
    this.outputService.executingBlocks.subscribe(
      (res) => { this.executingBlocks = res; },
    );
  }

  removeBlock(): void {
    this.blockService.removeBlock(this.block.blockId);
  }
}
