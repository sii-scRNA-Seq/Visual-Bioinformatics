import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { MatSnackBar, MatSnackBarModule } from '@angular/material/snack-bar';
import { TestBed } from '@angular/core/testing';
import { first, firstValueFrom } from 'rxjs';

import { BlockService } from './block.service';
import { CurrentDatasetService } from './current-dataset.service';
import { DatasetInfoService } from './dataset-info.service';
import { MockCurrentDatasetService } from './mock-current-dataset.service';
import { MockDatasetInfoService } from './mock-dataset-info.service';
import { MockOutputService } from './mock-output.service';
import { OutputService } from './output.service';

describe('BlockService', () => {
  let service: BlockService;
  let snackBar: MatSnackBar;
  let currentDatasetService: CurrentDatasetService;
  let datasetInfoService: DatasetInfoService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [
        BrowserAnimationsModule,
        MatSnackBarModule
      ],
      providers: [
        { provide: CurrentDatasetService, useClass: MockCurrentDatasetService },
        { provide: DatasetInfoService, useClass: MockDatasetInfoService },
        { provide: OutputService, useClass: MockOutputService }
      ]
    });
    service = TestBed.inject(BlockService);
    snackBar = TestBed.inject(MatSnackBar);
    currentDatasetService = TestBed.inject(CurrentDatasetService);
    datasetInfoService = TestBed.inject(DatasetInfoService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('BlocksOnCanvas', () => {
    it('should initially be empty', async () => {
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
    });

    it('should have the correct options on Load Data block before loading DatasetInfo', async () => {
      let blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas.length).toBe(0);
      service.addBlock('loaddata');
      blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas[0].blockId).toBe('loaddata');
      expect(blocksOnCanvas[0].parameters[0].type).toBe('SelectParameter');
      expect(blocksOnCanvas[0].parameters[0].value).toBe('');
      expect(blocksOnCanvas[0].parameters[0].options).toEqual([]);
    });

    it('should have the correct options on Load Data block after loading DatasetInfo', async () => {
      datasetInfoService.setDatasetInfo();
      let blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas.length).toBe(0);
      service.addBlock('loaddata');
      blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas[0].blockId).toBe('loaddata');
      expect(blocksOnCanvas[0].parameters[0].type).toBe('SelectParameter');
      expect(blocksOnCanvas[0].parameters[0].value).toBe('');
      expect(blocksOnCanvas[0].parameters[0].options).toEqual([
        {key: 'option1', text: 'Option 1'},
        {key: 'option2', text: 'Option 2'}
      ]);
    });

    it('should have the correct options on Load Data block after loading DatasetInfo and setting CurrentDataset', async () => {
      datasetInfoService.setDatasetInfo();
      currentDatasetService.setCurrentDataset('option1');
      let blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas.length).toBe(0);
      service.addBlock('loaddata');
      blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas[0].blockId).toBe('loaddata');
      expect(blocksOnCanvas[0].parameters[0].type).toBe('SelectParameter');
      expect(blocksOnCanvas[0].parameters[0].value).toBe('option1');
      expect(blocksOnCanvas[0].parameters[0].options).toEqual([
        {key: 'option1', text: 'Option 1'},
        {key: 'option2', text: 'Option 2'}
      ]);
    });

    it('should have the correct options on QC Filtering block before setting CurrentDataset', async () => {
      let blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas.length).toBe(0);
      service.addBlock('loaddata');
      service.addBlock('qcfiltering');
      blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas[1].blockId).toBe('qcfiltering');
      expect(blocksOnCanvas[1].parameters[0].type).toBe('SelectParameter');
      expect(blocksOnCanvas[1].parameters[0].value).toBe('');
      expect(blocksOnCanvas[1].parameters[0].options).toEqual([]);
    });

    it('should have the correct options on QC Filtering block after setting CurrentDataset', async () => {
      currentDatasetService.setCurrentDataset('option1');
      let blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas.length).toBe(0);
      service.addBlock('loaddata');
      service.addBlock('qcfiltering');
      blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas[1].blockId).toBe('qcfiltering');
      expect(blocksOnCanvas[1].parameters[0].type).toBe('SelectParameter');
      expect(blocksOnCanvas[1].parameters[0].value).toBe('sample1');
      expect(blocksOnCanvas[1].parameters[0].options).toEqual([
        {key: 'sample1', text: 'sample1'},
        {key: 'sample2', text: 'sample2'}
      ]);
    });

    it('should have the correct options on Integration block before setting CurrentDataset', async () => {
      let blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas.length).toBe(0);
      service.addBlock('loaddata');
      service.addBlock('variablegenes');
      service.addBlock('pca');
      service.addBlock('integration');
      blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas[3].blockId).toBe('integration');
      expect(blocksOnCanvas[3].parameters[0].type).toBe('SelectParameter');
      expect(blocksOnCanvas[3].parameters[0].value).toBe('');
      expect(blocksOnCanvas[3].parameters[0].options).toEqual([]);
    });

    it('should have the correct options on Integration block after setting CurrentDataset', async () => {
      currentDatasetService.setCurrentDataset('option1');
      let blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas.length).toBe(0);
      service.addBlock('loaddata');
      service.addBlock('variablegenes');
      service.addBlock('pca');
      service.addBlock('integration');
      blocksOnCanvas = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocksOnCanvas[3].blockId).toBe('integration');
      expect(blocksOnCanvas[3].parameters[0].type).toBe('SelectParameter');
      expect(blocksOnCanvas[3].parameters[0].value).toBe('ob1');
      expect(blocksOnCanvas[3].parameters[0].options).toEqual([
        {key: 'ob1', text: 'ob1'},
        {key: 'ob2', text: 'ob2'}
      ]);
    });
  });

  describe('addBlock', () => {
    beforeEach(() => {
      datasetInfoService.setDatasetInfo();
    });

    it('should add the given block when the ordering is valid - standard order', async () => {
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
      service.addBlock('loaddata');
      service.addBlock('basicfiltering');
      service.addBlock('qcplots');
      service.addBlock('qcfiltering');
      service.addBlock('variablegenes');
      service.addBlock('pca');
      service.addBlock('integration');
      service.addBlock('runumap');
      const result = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(result.length).toBe(8);
      expect(result[0].blockId).toBe('loaddata');
      expect(result[1].blockId).toBe('basicfiltering');
      expect(result[2].blockId).toBe('qcplots');
      expect(result[3].blockId).toBe('qcfiltering');
      expect(result[4].blockId).toBe('variablegenes');
      expect(result[5].blockId).toBe('pca');
      expect(result[6].blockId).toBe('integration');
      expect(result[7].blockId).toBe('runumap');
    });

    it('should add the given block when the ordering is valid - alternative order', async () => {
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
      service.addBlock('loaddata');
      service.addBlock('qcfiltering');
      service.addBlock('qcfiltering');
      service.addBlock('qcplots');
      service.addBlock('qcplots');
      service.addBlock('basicfiltering');
      service.addBlock('basicfiltering');
      service.addBlock('variablegenes');
      service.addBlock('variablegenes');
      service.addBlock('pca');
      service.addBlock('pca');
      service.addBlock('runumap');
      service.addBlock('runumap');
      service.addBlock('integration');
      service.addBlock('integration');
      service.addBlock('runumap');
      service.addBlock('runumap');
      const result = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(result.length).toBe(17);
      expect(result[0].blockId).toBe('loaddata');
      expect(result[1].blockId).toBe('qcfiltering');
      expect(result[2].blockId).toBe('qcfiltering');
      expect(result[3].blockId).toBe('qcplots');
      expect(result[4].blockId).toBe('qcplots');
      expect(result[5].blockId).toBe('basicfiltering');
      expect(result[6].blockId).toBe('basicfiltering');
      expect(result[7].blockId).toBe('variablegenes');
      expect(result[8].blockId).toBe('variablegenes');
      expect(result[9].blockId).toBe('pca');
      expect(result[10].blockId).toBe('pca');
      expect(result[11].blockId).toBe('runumap');
      expect(result[12].blockId).toBe('runumap');
      expect(result[13].blockId).toBe('integration');
      expect(result[14].blockId).toBe('integration');
      expect(result[15].blockId).toBe('runumap');
      expect(result[16].blockId).toBe('runumap');
    });

    it('should open snack bar when ordering is not valid - repeated load data', async() => {
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
      const spy = spyOn(snackBar, 'open');
      service.addBlock('loaddata');
      service.addBlock('loaddata');
      expect(spy).toHaveBeenCalledOnceWith('Load Data block cannot be added.', 'Close', { duration: 5000 });
    });
    
    it('should open snack bar when ordering is not valid - basic filtering before load data', async() => {
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
      const spy = spyOn(snackBar, 'open');
      service.addBlock('basicfiltering');
      expect(spy).toHaveBeenCalledOnceWith('Basic Filtering block cannot be added.', 'Close', { duration: 5000 });
    });

    it('should open snack bar when ordering is not valid - pca before variable genes', async() => {
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
      const spy = spyOn(snackBar, 'open');
      service.addBlock('loaddata');
      service.addBlock('pca');
      expect(spy).toHaveBeenCalledOnceWith('Principal Component Analysis block cannot be added.', 'Close', { duration: 5000 });
    });

    it('should open snack bar when ordering is not valid - integration before pca', async() => {
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
      const spy = spyOn(snackBar, 'open');
      service.addBlock('loaddata');
      service.addBlock('variablegenes');
      service.addBlock('integration');
      expect(spy).toHaveBeenCalledOnceWith('Integration block cannot be added.', 'Close', { duration: 5000 });
    });

    it('should open snack bar when ordering is not valid - run umap before pca', async() => {
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
      const spy = spyOn(snackBar, 'open');
      service.addBlock('loaddata');
      service.addBlock('variablegenes');
      service.addBlock('runumap');
      expect(spy).toHaveBeenCalledOnceWith('Run UMAP block cannot be added.', 'Close', { duration: 5000 });
    });
  });

  describe('removeBlock', () => {
    beforeEach(() => {
      datasetInfoService.setDatasetInfo();
    });

    it('should remove the only block when called for the only block', async () => {
      service.addBlock('loaddata');
      let blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(1);
      service.removeBlock(blocks[0].blockUUID);
      blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
    });

    it('should remove the final block when called for the final block', async () => {
      service.addBlock('loaddata');
      service.addBlock('basicfiltering');
      let blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(2);
      service.removeBlock(blocks[1].blockUUID);
      blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(1);
      expect(blocks[0].blockId).toBe('loaddata');
    });

    it('should remove all blocks after the given block when called for a non-final block', async () => {
      service.addBlock('loaddata');
      service.addBlock('basicfiltering');
      service.addBlock('basicfiltering');
      let blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(3);
      service.removeBlock(blocks[1].blockUUID);
      blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(1);
      expect(blocks[0].blockId).toBe('loaddata');
    });

    it('should remove only the given version when multiple blocks share a blockId', async () => {
      service.addBlock('loaddata');
      service.addBlock('basicfiltering');
      service.addBlock('basicfiltering');
      let blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(3);
      service.removeBlock(blocks[2].blockUUID);
      blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(2);
      expect(blocks[0].blockId).toBe('loaddata');
      expect(blocks[1].blockId).toBe('basicfiltering');
    });

    it('should have no effect when no blocks have the given blockUUID', async () => {
      service.addBlock('loaddata');
      service.addBlock('basicfiltering');
      service.addBlock('basicfiltering');
      let blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(3);
      service.removeBlock('ImpossibleUUID');
      blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(3);
    });
  });

  describe('executeBlocks', () => {
    it('should result in a call of outputService.resetOutputs', async () => {
      const outputService: OutputService = TestBed.inject(OutputService);
      spyOn(outputService, 'resetOutputs');
      service.executeBlocks();
      expect(outputService.resetOutputs).toHaveBeenCalledTimes(1);
    });

    it('should result in a call of outputService.executeBlocks', async () => {
      const outputService: OutputService = TestBed.inject(OutputService);
      spyOn(outputService, 'executeBlocks');
      service.executeBlocks();
      expect(outputService.executeBlocks).toHaveBeenCalledTimes(1);
    });
  });
});
