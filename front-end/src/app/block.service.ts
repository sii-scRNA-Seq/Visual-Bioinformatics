import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';
import { MatSnackBar } from '@angular/material/snack-bar';

import { Block, BlockId } from './block.interface';
import { OutputService } from './output.service';

@Injectable({
  providedIn: 'root',
})

export class BlockService {
  private readonly blocksOnCanvas$: BehaviorSubject<Block[]> = new BehaviorSubject<Block[]> ([]);
  readonly blocksOnCanvas: Observable<Block[]> = this.blocksOnCanvas$.asObservable();

  constructor(private outputService: OutputService, private snackBar: MatSnackBar) { }

  addBlock(id: BlockId): void {
    const blockList = this.blocksOnCanvas$.getValue();
    switch (id) {
      case 'loaddata': {
        if (blockList.length == 0) {
          this.blocksOnCanvas$.next([{
            blockId: 'loaddata',
            title: 'Load Data',
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering','variablegenes'],
            parameters: [],
          }]);
        }
        else {
          this.snackBar.open('Load Data block cannot be added', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'basicfiltering': {
        if (blockList[blockList.length-1]?.possibleChildBlocks.indexOf('basicfiltering') > -1) {
          blockList.push({
            blockId: 'basicfiltering',
            title: 'Basic Filtering',
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering','variablegenes'],
            parameters: [
              {key: 'min_genes', text: 'Minimum Genes Per Cell', value: 200},
              {key: 'min_cells', text: 'Minimum Cells Per Gene', value: 3}
            ],
          });
          this.blocksOnCanvas$.next(blockList);
        }
        else {
          this.snackBar.open('Basic Filtering block cannot be added', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'qcplots': {
        if (blockList[blockList.length-1]?.possibleChildBlocks.indexOf('qcplots') > -1) {
          blockList.push({
            blockId: 'qcplots',
            title: 'Quality Control Plots',
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering','variablegenes'],
            parameters: [],
          });
          this.blocksOnCanvas$.next(blockList);
        }
        else {
          this.snackBar.open('Quality Control Plots block cannot be added', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'qcfiltering': {
        if (blockList[blockList.length-1]?.possibleChildBlocks.indexOf('qcfiltering') > -1) {
          blockList.push({
            blockId: 'qcfiltering',
            title: 'Quality Control Filtering',
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering','variablegenes'],
            parameters: [
              {key: 'n_genes_by_counts', text: 'Maximum Total Gene Count', value: 2500},
              {key: 'pct_counts_mt', text: 'Maximum % Mitochondrial Genes', value: 5}
            ],
          });
          this.blocksOnCanvas$.next(blockList);
        }
        else {
          this.snackBar.open('Quality Control Filtering block cannot be added', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'variablegenes': {
        if (blockList[blockList.length-1]?.possibleChildBlocks.indexOf('variablegenes') > -1) {
          blockList.push({
            blockId: 'variablegenes',
            title: 'Identify Highly Variable Genes',
            possibleChildBlocks: ['variablegenes', 'pca'],
            parameters: [
              {key: 'min_mean', text: 'Minimum Mean', value: 0.0125},
              {key: 'max_mean', text: 'Maximum Mean', value: 3},
              {key: 'min_disp', text: 'Minimum Dispersion', value: 0.5}
            ],
          });
          this.blocksOnCanvas$.next(blockList);
        }
        else {
          this.snackBar.open('Identify Highly Variable Genes block cannot be added', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'pca': {
        if (blockList[blockList.length-1]?.possibleChildBlocks.indexOf('pca') > -1) {
          blockList.push({
            blockId: 'pca',
            title: 'Principle Component Analysis',
            possibleChildBlocks: ['pca'],
            parameters: [],
          });
          this.blocksOnCanvas$.next(blockList);
        }
        else {
          this.snackBar.open('Principle Component Analysis block cannot be added', 'Close', { duration: 5000 });
        }
        break;
      }
    }
  }

  removeBlock(id: BlockId): void {
    const newBlockList: Block[] = [];
    for (let i = 0; i < this.blocksOnCanvas$.getValue().length; i++) {
      if (this.blocksOnCanvas$.getValue()[i].blockId == id) {
        this.blocksOnCanvas$.next(newBlockList);
        break;
      }
      else {
        newBlockList.push(this.blocksOnCanvas$.getValue()[i]);
      }
    }
  }

  async executeBlocks(): Promise<void> {
    this.outputService.resetOutputs();
    for(let i=0; i < this.blocksOnCanvas$.getValue().length; i++) {
      await this.outputService.executeBlock(this.blocksOnCanvas$.getValue()[i]);
    }
  }
}
