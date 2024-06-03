import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';
import { MatSnackBar } from '@angular/material/snack-bar';

import { Block, BlockId } from './block.interface';
import { BlockServiceInterface } from './block.service.interface';
import { OutputService } from './output.service';

@Injectable({
  providedIn: 'root',
})
export class BlockService implements BlockServiceInterface {
  private readonly blocksOnCanvas$: BehaviorSubject<Block[]> = new BehaviorSubject<Block[]> ([]);
  readonly blocksOnCanvas: Observable<Block[]> = this.blocksOnCanvas$.asObservable();

  constructor(private outputService: OutputService, private snackBar: MatSnackBar) { }

  addBlock(id: BlockId): void {
    const blockList = this.blocksOnCanvas$.getValue();
    const lastBlock = blockList[blockList.length-1];
    switch (id) {
      case 'loaddata': {
        if (blockList.length == 0) {
          this.blocksOnCanvas$.next([{
            blockId: 'loaddata',
            title: 'Load Data',
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering','variablegenes'],
            parameters: [
              {key: 'dataset', text: 'Dataset', options: ['A', 'B', 'C'], value: 'A'},
            ],
          }]);
        }
        else {
          this.snackBar.open('Load Data block cannot be added.', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'basicfiltering': {
        if (lastBlock?.possibleChildBlocks.indexOf('basicfiltering') > -1) {
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
          this.snackBar.open('Basic Filtering block cannot be added.', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'qcplots': {
        if (lastBlock?.possibleChildBlocks.indexOf('qcplots') > -1) {
          blockList.push({
            blockId: 'qcplots',
            title: 'Quality Control Plots',
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering','variablegenes'],
            parameters: [],
          });
          this.blocksOnCanvas$.next(blockList);
        }
        else {
          this.snackBar.open('Quality Control Plots block cannot be added.', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'qcfiltering': {
        if (lastBlock?.possibleChildBlocks.indexOf('qcfiltering') > -1) {
          blockList.push({
            blockId: 'qcfiltering',
            title: 'Quality Control Filtering',
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering','variablegenes'],
            parameters: [
              {key: 'min_n_genes_by_counts', text: 'Minimum Genes Per Cell', value: 200},
              {key: 'max_n_genes_by_counts', text: 'Maximum Gene Per Cell', value: 2500},
              {key: 'pct_counts_mt', text: 'Maximum % Mitochondrial Genes', value: 5}
            ],
          });
          this.blocksOnCanvas$.next(blockList);
        }
        else {
          this.snackBar.open('Quality Control Filtering block cannot be added.', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'variablegenes': {
        if (lastBlock?.possibleChildBlocks.indexOf('variablegenes') > -1) {
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
          this.snackBar.open('Identify Highly Variable Genes block cannot be added.', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'pca': {
        if (lastBlock?.possibleChildBlocks.indexOf('pca') > -1) {
          blockList.push({
            blockId: 'pca',
            title: 'Principal Component Analysis',
            possibleChildBlocks: ['pca', 'runumap'],
            parameters: [],
          });
          this.blocksOnCanvas$.next(blockList);
        }
        else {
          this.snackBar.open('Principal Component Analysis block cannot be added.', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'runumap': {
        if (lastBlock?.possibleChildBlocks.indexOf('runumap') > -1) {
          blockList.push({
            blockId: 'runumap',
            title: 'Run UMAP',
            possibleChildBlocks: ['runumap'],
            parameters: [
              {key: 'n_neighbors', text: 'Number of Neighbours', value: 10},
              {key: 'n_pcs', text: 'Number of Principal Components', value: 40},
            ],
          });
          this.blocksOnCanvas$.next(blockList);
        }
        else {
          this.snackBar.open('Run UMAP block cannot be added.', 'Close', { duration: 5000 });
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

  executeBlocks(): void {
    this.outputService.resetOutputs();
    this.outputService.executeBlocks(this.blocksOnCanvas$.getValue());
  }
}
