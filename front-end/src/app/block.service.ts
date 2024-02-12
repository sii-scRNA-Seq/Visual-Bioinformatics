import { BehaviorSubject, from, Observable } from 'rxjs';
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
    switch (id) {
      case 'loaddata': {
        if (this.blocksOnCanvas$.getValue().length == 0) {
          this.blocksOnCanvas$.next([{
            blockId: 'loaddata',
            title: 'Load Data',
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering'],
            parameters: [],
            onRun: () => from(''),
          }]);
        }
        else {
          this.snackBar.open('Load Data block cannot be added', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'basicfiltering': {
        if (this.blocksOnCanvas$.getValue()[this.blocksOnCanvas$.getValue().length-1]?.possibleChildBlocks.indexOf('basicfiltering') > -1) {
          const temp = this.blocksOnCanvas$.getValue();
          temp.push({
            blockId: 'basicfiltering',
            title: 'Basic Filtering',
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering'],
            parameters: [
              {key: 'min_genes', text: 'Minimum Genes Per Cell', value: 200},
              {key: 'min_cells', text: 'Minimum Cells Per Gene', value: 3}
            ],
            onRun: () => from(''),
          });
          this.blocksOnCanvas$.next(temp);
        }
        else {
          this.snackBar.open('Basic Filtering block cannot be added', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'qcplots': {
        if (this.blocksOnCanvas$.getValue()[this.blocksOnCanvas$.getValue().length-1]?.possibleChildBlocks.indexOf('qcplots') > -1) {
          const temp = this.blocksOnCanvas$.getValue();
          temp.push({
            blockId: 'qcplots',
            title: 'Quality Control Plots',
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering'],
            parameters: [],
            onRun: () => from(''),
          });
          this.blocksOnCanvas$.next(temp);
        }
        else {
          this.snackBar.open('Quality Control Plots block cannot be added', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'qcfiltering': {
        if (this.blocksOnCanvas$.getValue()[this.blocksOnCanvas$.getValue().length-1]?.possibleChildBlocks.indexOf('qcfiltering') > -1) {
          const temp = this.blocksOnCanvas$.getValue();
          temp.push({
            blockId: 'qcfiltering',
            title: 'Quality Control Filtering',
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering'],
            parameters: [
              {key: 'n_genes_by_counts', text: 'Maximum Total Gene Count', value: 2500},
              {key: 'pct_counts_mt', text: 'Maximum % Mitochondrial Genes', value: 5}
            ],
            onRun: () => from(''),
          });
          this.blocksOnCanvas$.next(temp);
        }
        else {
          this.snackBar.open('Quality Control Filtering block cannot be added', 'Close', { duration: 5000 });
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
