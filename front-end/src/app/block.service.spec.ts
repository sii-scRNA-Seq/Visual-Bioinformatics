import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { first, firstValueFrom } from 'rxjs';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatSnackBar, MatSnackBarModule } from '@angular/material/snack-bar';
import { TestBed } from '@angular/core/testing';

import { BlockService } from './block.service';
import { DatasetInfoService } from './dataset-info.service';
import { MockDatasetInfoService } from './mock-dataset-info.service';
import { MockOutputService } from './mock-output.service';
import { OutputService } from './output.service';

describe('BlockService', () => {
  let service: BlockService;
  let snackBar: MatSnackBar;
  let datasetInfoService: DatasetInfoService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [
        BrowserAnimationsModule,
        HttpClientTestingModule,
        MatSnackBarModule,
      ],
      providers: [
        { provide: OutputService, useClass: MockOutputService },
        { provide: DatasetInfoService, useClass: MockDatasetInfoService },
      ],
    });
    service = TestBed.inject(BlockService);
    snackBar = TestBed.inject(MatSnackBar);
    datasetInfoService = TestBed.inject(DatasetInfoService);
    datasetInfoService.setDatasetInfo();
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('BlocksOnCanvas', () => {
    it('should initially be empty', async () => {
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
    });
  });

  describe('addBlock', () => {
    it('should add the given block when the ordering is valid - standard order', async () => {
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
      service.addBlock('loaddata');
      service.addBlock('basicfiltering');
      service.addBlock('qcplots');
      service.addBlock('qcfiltering');
      service.addBlock('variablegenes');
      service.addBlock('pca');
      service.addBlock('runumap');
      const result = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(result.length).toBe(7);
      expect(result[0].blockId).toBe('loaddata');
      expect(result[1].blockId).toBe('basicfiltering');
      expect(result[2].blockId).toBe('qcplots');
      expect(result[3].blockId).toBe('qcfiltering');
      expect(result[4].blockId).toBe('variablegenes');
      expect(result[5].blockId).toBe('pca');
      expect(result[6].blockId).toBe('runumap');
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
      const result = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(result.length).toBe(13);
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
    it('should remove the only block when called for only one block', async () => {
      service.addBlock('loaddata');
      const blocks1 = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks1.length).toBe(1);
      service.removeBlock('loaddata');
      const results1 = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(results1.length).toBe(0);
    });

    it('should remove the final block when called for the final block', async () => {
      service.addBlock('loaddata');
      service.addBlock('basicfiltering');
      const blocks2 = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks2.length).toBe(2);
      service.removeBlock('basicfiltering');
      const results2 = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(results2.length).toBe(1);
      expect(results2[0].blockId).toBe('loaddata');
    });

    it('should remove all blocks after the given block when called for a non-final block', async () => {
      service.addBlock('loaddata');
      service.addBlock('basicfiltering');
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(2);
      service.removeBlock('loaddata');
      const results = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(results.length).toBe(0);
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
