import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';
import { MatSnackBar } from '@angular/material/snack-bar';

import { Block, BlockId, Option } from './block.interface';
import { BlockServiceInterface } from './block.service.interface';
import { DatasetInfo } from './dataset-info';
import { DatasetInfoService } from './dataset-info.service';
import { OutputService } from './output.service';

@Injectable({
  providedIn: 'root',
})
export class BlockService implements BlockServiceInterface {
  private readonly blocksOnCanvas$: BehaviorSubject<Block[]> = new BehaviorSubject<Block[]> ([]);
  readonly blocksOnCanvas: Observable<Block[]> = this.blocksOnCanvas$.asObservable();

  private datasetInfo: DatasetInfo[] = [];

  private currentDataset: string | number = '';

  constructor(private outputService: OutputService, private snackBar: MatSnackBar, private datasetInfoService: DatasetInfoService) { 
    this.datasetInfoService.datasetInfo.subscribe(
      (res) => {
        this.datasetInfo = res;
        if (this.datasetInfo.length > 0) {
          this.currentDataset = this.datasetInfo[0].key;
        }
      },
    );
    this.blocksOnCanvas.subscribe(
      (res) => {
        if (res.length > 0) {
          this.currentDataset = res[0].parameters[0].value || '';
        }
      }
    );
  }

  addBlock(id: BlockId): void {
    const blockList = this.blocksOnCanvas$.getValue();
    const lastBlock = blockList[blockList.length-1];
    switch (id) {
      case 'loaddata': {
        if (blockList.length == 0) {
          const options: Option[] = [];
          this.datasetInfo.forEach( dataset => {
            options.push({key: dataset.key, text: dataset.title});
          });
          this.blocksOnCanvas$.next([{
            blockId: 'loaddata',
            title: 'Load Data',
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering','variablegenes'],
            parameters: [
              {type: 'SelectParameter', key: 'dataset', text: 'Dataset', value: this.currentDataset, options: options},
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
              {type: 'InputParameter', key: 'min_genes', text: 'Minimum Genes Per Cell', value: 200},
              {type: 'InputParameter', key: 'min_cells', text: 'Minimum Cells Per Gene', value: 3}
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
              {type: 'InputParameter', key: 'min_n_genes_by_counts', text: 'Minimum Genes Per Cell', value: 200},
              {type: 'InputParameter', key: 'max_n_genes_by_counts', text: 'Maximum Gene Per Cell', value: 2500},
              {type: 'InputParameter', key: 'pct_counts_mt', text: 'Maximum % Mitochondrial Genes', value: 5}
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
              {type: 'InputParameter', key: 'min_mean', text: 'Minimum Mean', value: 0.0125},
              {type: 'InputParameter', key: 'max_mean', text: 'Maximum Mean', value: 3},
              {type: 'InputParameter', key: 'min_disp', text: 'Minimum Dispersion', value: 0.5}
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
            possibleChildBlocks: ['pca', 'integration', 'runumap'],
            parameters: [],
          });
          this.blocksOnCanvas$.next(blockList);
        }
        else {
          this.snackBar.open('Principal Component Analysis block cannot be added.', 'Close', { duration: 5000 });
        }
        break;
      }
      case 'integration': {
        if (lastBlock?.possibleChildBlocks.indexOf('integration') > -1) {
          let value: string = '';
          const options: Option[] = [];
          this.datasetInfo.forEach( dataset => {
            if (dataset.key == this.currentDataset && dataset.integration_obs.length > 0) {
              value = dataset.integration_obs[0];
              options.length = 0;
              dataset.integration_obs.forEach( observation => {
                options.push({key: observation, text: observation});
              });
            }
          });

          blockList.push({
            blockId: 'integration',
            title: 'Integration',
            possibleChildBlocks: ['integration', 'runumap'],
            parameters: [
              {type: 'SelectParameter', key: 'observation', text: 'Observation', value: value, options: options},
            ],
          });
          this.blocksOnCanvas$.next(blockList);
        }
        else {
          this.snackBar.open('Integration block cannot be added.', 'Close', { duration: 5000 });
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
              {type: 'InputParameter', key: 'n_neighbors', text: 'Number of Neighbours', value: 10},
              {type: 'InputParameter', key: 'n_pcs', text: 'Number of Principal Components', value: 40},
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
