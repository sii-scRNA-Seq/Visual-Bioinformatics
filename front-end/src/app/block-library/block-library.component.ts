import { Component, OnInit } from '@angular/core';

import { Block, BlockId } from '../block.interface';
import { BlockService } from '../block.service';
import { OutputService } from '../output.service';

@Component({
  selector: 'app-block-library',
  templateUrl: './block-library.component.html',
  styleUrls: ['./block-library.component.css'],
})

export class BlockLibraryComponent implements OnInit {
  blockList: Block[] = [];

  executingBlocks: boolean = false;

  blockLibrary: BlockStatus[] = [
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

  ngOnInit() {
    this.updateDisabledBlocks();
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

  getTooltipMessage(): string {
    if (this.executingBlocks) {
      return 'Blocks cannot be added while blocks are being executed.';
    } else if (this.blockList.length == 0) {
      return 'This block cannot be added to an empty canvas.';
    } else if (this.blockList[this.blockList.length - 1].title.match('^[AEIOU].*')) {
      return 'This block cannot be immediately below an ' + this.blockList[this.blockList.length - 1].title + ' block.';
    } else {
      return 'This block cannot be immediately below a ' + this.blockList[this.blockList.length - 1].title + ' block.';
    }
  }
}

interface BlockStatus {
  blockId: BlockId
  title: string
  disabled: boolean
}
