import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';
import { MatSnackBar } from '@angular/material/snack-bar';
import { v4 as uuidv4 } from 'uuid';

import {Block, BlockId, BlockIdToTitleMap, Option} from './block.interface';
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

  private currentDataset: DatasetInfo = {
    key: '',
    title: '',
    integration_obs: []
  };

  constructor(private outputService: OutputService, private snackBar: MatSnackBar, private datasetInfoService: DatasetInfoService) {
    this.datasetInfoService.datasetInfo.subscribe(
      (res) => {
        this.datasetInfo = res;
        if (this.datasetInfo.length > 0) {
          this.currentDataset = this.datasetInfo[0];
        }
      },
    );
    this.blocksOnCanvas.subscribe(
      (res) => {
        if (res.length > 0) {
          this.currentDataset = this.datasetInfo.find(dataset => dataset.key == res[0].parameters[0].value) || this.datasetInfo[0] || this.currentDataset;
        }
      }
    );
  }

  onMatSelectValueChanges(block: Block): void {
    if (block.blockId == 'loaddata') {
      this.currentDataset = this.datasetInfo.find(dataset => dataset.key == block.parameters[0].value) || this.datasetInfo[0] || this.currentDataset;
      const blockList = this.blocksOnCanvas$.getValue();
      for (let i = 0; i < blockList.length; i++) {
        if (blockList[i].blockId == 'integration') {
          let value: string = '';
          const options: Option[] = [];
          if (this.currentDataset.integration_obs.length > 0) {
            value = this.currentDataset.integration_obs[0];
            this.currentDataset.integration_obs.forEach( observation => {
              options.push({key: observation, text: observation});
            });
          }
          blockList[i].parameters[0].options = options;
          blockList[i].parameters[0].value = value;          
        }
      }
      console.log(this.currentDataset.key);
    } else {
      console.log("ERROR");
    }
  }

  addBlock(id: BlockId): void {
    const blockList = this.blocksOnCanvas$.getValue();
    const lastBlock: Block = blockList[blockList.length-1];
    switch (id) {
      case 'loaddata': {
        if (blockList.length == 0) {
          const options: Option[] = [];
          this.datasetInfo.forEach( dataset => {
            options.push({key: dataset.key, text: dataset.title});
          });
          this.blocksOnCanvas$.next([{
            blockId: 'loaddata',
            blockUUID: uuidv4(),
            title: BlockIdToTitleMap.loaddata,
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering','variablegenes'],
            parameters: [
              {type: 'SelectParameter', key: 'dataset', text: 'Dataset', value: this.currentDataset.key, options: options},
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
            blockUUID: uuidv4(),
            title: BlockIdToTitleMap.basicfiltering,
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
            blockUUID: uuidv4(),
            title: BlockIdToTitleMap.qcplots,
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
            blockUUID: uuidv4(),
            title: BlockIdToTitleMap.qcfiltering,
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
            blockUUID: uuidv4(),
            title: BlockIdToTitleMap.variablegenes,
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
            blockUUID: uuidv4(),
            title: BlockIdToTitleMap.pca,
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
          if (this.currentDataset.integration_obs.length > 0) {
            value = this.currentDataset.integration_obs[0];
            this.currentDataset.integration_obs.forEach( observation => {
              options.push({key: observation, text: observation});
            });
          }
          blockList.push({
            blockId: 'integration',
            blockUUID: uuidv4(),
            title: BlockIdToTitleMap.integration,
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
            blockUUID: uuidv4(),
            title: BlockIdToTitleMap.runumap,
            possibleChildBlocks: ['integration', 'runumap'],
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

  removeBlock(uuid: string): void {
    const newBlockList: Block[] = [];
    for (let i = 0; i < this.blocksOnCanvas$.getValue().length; i++) {
      if (this.blocksOnCanvas$.getValue()[i].blockUUID == uuid) {
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
