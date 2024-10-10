import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';
import { MatSnackBar } from '@angular/material/snack-bar';
import { v4 as uuidv4 } from 'uuid';

import {Block, BlockId, BlockIdToTitleMap, Option} from './block.interface';
import { BlockServiceInterface } from './block.service.interface';
import { CurrentDatasetService } from './current-dataset.service';
import { DatasetInfo } from './dataset-info';
import { DatasetInfoService } from './dataset-info.service';
import { OutputService } from './output.service';

@Injectable({
  providedIn: 'root',
})
export class BlockService implements BlockServiceInterface {
  private readonly blocksOnCanvas$: BehaviorSubject<Block[]> = new BehaviorSubject<Block[]> ([]);
  readonly blocksOnCanvas: Observable<Block[]> = this.blocksOnCanvas$.asObservable();

  private currentDataset: DatasetInfo = {
    key: '',
    title: '',
    samples: [],
    integration_obs: []
  };

  private datasetInfo: DatasetInfo[] = [];

  constructor(private outputService: OutputService, private snackBar: MatSnackBar, private datasetInfoService: DatasetInfoService, private currentDatasetService: CurrentDatasetService) {
    this.currentDatasetService.currentDataset.subscribe(
      (res) => {
        this.currentDataset = res;
        this.onCurrentDatasetChanges();
      },
    );
    this.datasetInfoService.datasetInfo.subscribe(
      (res) => {
        this.datasetInfo = res;
      },
    );
  }

  /**
   * Updates the value and options for the sample parameter on all 'qcfiltering' blocks to match the current dataset.
   * Updates the value and options for the observation parameter on all 'integration' blocks to match the current dataset.
   */
  onCurrentDatasetChanges(): void {
    const blockList = this.blocksOnCanvas$.getValue();
    for (let i = 0; i < blockList.length; i++) {
      if (blockList[i].blockId == 'qcfiltering') {
        let sampleValue: string = '';
        const sampleOptions: Option[] = [];
        if (this.currentDataset.samples.length > 0) {
          sampleValue = this.currentDataset.samples[0];
          this.currentDataset.samples.forEach( sample => {
            sampleOptions.push({key: sample, text: sample});
          });
        }
        blockList[i].parameters[0].options = sampleOptions;
        blockList[i].parameters[0].value = sampleValue;
      } else if (blockList[i].blockId == 'integration') {
        let observationValue: string = '';
        const observationOptions: Option[] = [];
        if (this.currentDataset.integration_obs.length > 0) {
          observationValue = this.currentDataset.integration_obs[0];
          this.currentDataset.integration_obs.forEach( observation => {
            observationOptions.push({key: observation, text: observation});
          });
        }
        blockList[i].parameters[0].options = observationOptions;
        blockList[i].parameters[0].value = observationValue;
      }
    }
  }

  /**
   * Creates a block with all the required information and adds it to the Canvas.
   *
   * @param id - The BlockId of the block that should be added
   */
  addBlock(id: BlockId): void {
    const blockList = this.blocksOnCanvas$.getValue();
    const lastBlock: Block = blockList[blockList.length-1];
    switch (id) {
      case 'loaddata': {
        if (blockList.length == 0) {
          const datasetOptions: Option[] = [];
          this.datasetInfo.forEach( dataset => {
            datasetOptions.push({key: dataset.key, text: dataset.title});
          });
          this.blocksOnCanvas$.next([{
            blockId: 'loaddata',
            blockUUID: uuidv4(),
            title: BlockIdToTitleMap.loaddata,
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering','variablegenes'],
            parameters: [
              {type: 'SelectParameter', key: 'dataset', text: 'Dataset', value: this.currentDataset.key, options: datasetOptions},
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
          const sampleValue: string = this.currentDataset.samples[0] || '';
          const sampleOptions: Option[] = [];
          this.currentDataset.samples.forEach( sample => {
            sampleOptions.push({key: sample, text: sample});
          });
          blockList.push({
            blockId: 'qcfiltering',
            blockUUID: uuidv4(),
            title: BlockIdToTitleMap.qcfiltering,
            possibleChildBlocks: ['basicfiltering','qcplots','qcfiltering','variablegenes'],
            parameters: [
              {type: 'SelectParameter', key: 'sample', text: 'Sample', value: sampleValue, options: sampleOptions},
              {type: 'InputParameter', key: 'min_n_genes_by_counts', text: 'Minimum Genes Per Cell', value: 200},
              {type: 'InputParameter', key: 'max_n_genes_by_counts', text: 'Maximum Genes Per Cell', value: 2500},
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
            possibleChildBlocks: ['pca', 'integration', 'runumap', 'plot_reddim'],
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
          const observationValue: string = this.currentDataset.integration_obs[0] || '';
          const observationOptions: Option[] = [];
          this.currentDataset.integration_obs.forEach( observation => {
            observationOptions.push({key: observation, text: observation});
          });
          blockList.push({
            blockId: 'integration',
            blockUUID: uuidv4(),
            title: BlockIdToTitleMap.integration,
            possibleChildBlocks: ['integration', 'runumap', 'plot_reddim'],
            parameters: [
              {type: 'SelectParameter', key: 'observation', text: 'Observation', value: observationValue, options: observationOptions},
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
            possibleChildBlocks: ['integration', 'runumap', 'plot_reddim'],
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
      case 'plot_reddim': {
        if (lastBlock?.possibleChildBlocks.indexOf('plot_reddim') > -1) {
          blockList.push({
            blockId: 'plot_reddim',
            blockUUID: uuidv4(),
            title: BlockIdToTitleMap.plot_reddim,
            possibleChildBlocks: ['integration', 'runumap', 'plot_reddim'],
            parameters: [
              {type: 'SelectParameter', key: 'reduction', text: 'Reduction', value: 'PCA', options: [
                {key: 'PCA', text: 'PCA'},
                {key: 'UMAP', text: 'UMAP'}
              ]},
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

  /**
   * Finds the first block with the given uuid and removes that block and all blocks below it from the Canvas.
   *
   * @param uuid - The uuid of the block that should be removed
   */
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

  /**
   * Calls the resetOutputs() method in the outputService to clear the Output Display.
   * Calls the executeBlocks() method in the outputService with all of the information of the blocks currently on the Canvas.
   */
  executeBlocks(): void {
    this.outputService.resetOutputs();
    this.outputService.executeBlocks(this.blocksOnCanvas$.getValue());
  }
}
