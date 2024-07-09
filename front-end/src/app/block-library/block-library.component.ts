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

  blockLibrary: LibraryInfo[] = [
    {blockId: 'loaddata', title: 'Load Data', disabled: true},
    {blockId: 'basicfiltering', title: 'Basic Filtering', disabled: true},
    {blockId: 'qcplots', title: 'Quality Control Plots', disabled: true},
    {blockId: 'qcfiltering', title: 'Quality Controls Filtering', disabled: true},
    {blockId: 'variablegenes', title: 'Identify Highly Variable Genes', disabled: true},
    {blockId: 'pca', title: 'Principal Component Analysis', disabled: true},
    {blockId: 'integration', title: 'Integration', disabled: true},
    {blockId: 'runumap', title: 'Run UMAP', disabled: true}
  ];

  constructor(private blockService: BlockService, private outputService: OutputService) {
    this.updateDisabledBlocks();
    this.blockService.blocksOnCanvas.subscribe(
      (res) => { 
        this.blockList = res;
        this.updateDisabledBlocks();
      },
    );
    this.outputService.executingBlocks.subscribe(
      (res) => {
        this.executingBlocks = res;
        this.updateDisabledBlocks();
      },
    );
  }

  addBlock(id: BlockId): void {
    this.blockService.addBlock(id);
  }

  updateDisabledBlocks(): void {
    this.blockLibrary[0].disabled = this.executingBlocks || this.blockList.length !== 0;
    for (let i = 1; i < this.blockLibrary.length; i++) {
      this.blockLibrary[i].disabled = this.executingBlocks || this.blockList.length === 0 || this.blockList[this.blockList.length - 1].possibleChildBlocks.indexOf(this.blockLibrary[i].blockId) < 0;
    }
  }
}

interface LibraryInfo {
  blockId: BlockId
  title: string
  disabled: boolean
}