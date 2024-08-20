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
    {blockId: 'loaddata', title: 'Load Data', disabled: true, helpText: 'This step allows us to read our data into SCAMPI ready for analysis!'},
    {blockId: 'basicfiltering', title: 'Basic Filtering', disabled: true, helpText: 'Remove cells and genes that are not expressing much'},
    {blockId: 'qcplots', title: 'Quality Control Plots', disabled: true, helpText: 'Plot some metrics to see if we have any outliers in our data'},
    {blockId: 'qcfiltering', title: 'Quality Controls Filtering', disabled: true, helpText: 'Apply a threshold to remove outliers that could be doublets or stressed/dying cells'},
    {blockId: 'variablegenes', title: 'Identify Highly Variable Genes', disabled: true, helpText: 'This selects a set of genes that explain the most variation in our data and explain most of the underlying biology'},
    {blockId: 'pca', title: 'Principal Component Analysis', disabled: true, helpText: 'PCA helps us group genes together in principle components that explain the biological differences in our data'},
    {blockId: 'integration', title: 'Integration', disabled: true, helpText: 'Unwanted sources of variation need to be corrected for, revealing differences that are driven by the biology and not batch effects'},
    {blockId: 'runumap', title: 'Run UMAP', disabled: true, helpText: 'This allows us to visualise our complex data in a simpler way'}
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

  /**
   * Calls the addBlock() method in the blockService with the blockId of the block to be added.
   * 
   * @param id - The blockId of the block to be added
   */
  addBlock(id: BlockId): void {
    this.blockService.addBlock(id);
  }

  /**
   * Updates the block library to show which blocks should be available/unavailable to be added at any given time.
   */
  updateDisabledBlocks(): void {
    this.blockLibrary[0].disabled = this.executingBlocks || this.blockList.length !== 0;
    for (let i = 1; i < this.blockLibrary.length; i++) {
      this.blockLibrary[i].disabled = this.executingBlocks || this.blockList.length === 0 || this.blockList[this.blockList.length - 1].possibleChildBlocks.indexOf(this.blockLibrary[i].blockId) < 0;
    }
  }

  /**
   * Chooses the text that should be displayed on tooltips at any given time.
   * 
   * @returns The text to display on the tooltip
   */
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
  helpText: string
}
