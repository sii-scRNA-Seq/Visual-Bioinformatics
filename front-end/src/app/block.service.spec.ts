import { first, firstValueFrom } from 'rxjs';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatSnackBarModule } from '@angular/material/snack-bar';
import { TestBed } from '@angular/core/testing';

import { BlockService } from './block.service';
import { OutputService } from './output.service';

describe('BlockService', () => {
  let service: BlockService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [
        HttpClientTestingModule,
        MatSnackBarModule,
      ],
    });
    service = TestBed.inject(BlockService);
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
    it('should add the given block when called', async () => {
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
      service.addBlock('loaddata');
      const results1 = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(results1.length).toBe(1);
      expect(results1[0].blockId).toBe('loaddata');
      service.addBlock('basicfiltering');
      const results2 = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(results2.length).toBe(2);
      expect(results2[1].blockId).toBe('basicfiltering');
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

    it('should result in no calls of outputService.executeBlock when no blocks are on the canvas', async () => {
      const outputService: OutputService = TestBed.inject(OutputService);
      spyOn(outputService, 'executeBlock');
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
      service.executeBlocks();
      expect(outputService.executeBlock).toHaveBeenCalledTimes(0);
    });

    it('should result in one call of outputService.executeBlock when one block is on the canvas', async () => {
      const outputService: OutputService = TestBed.inject(OutputService);
      spyOn(outputService, 'executeBlock');
      service.addBlock('loaddata');
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(1);
      service.executeBlocks();
      expect(outputService.executeBlock).toHaveBeenCalledTimes(1);
    });

    it('should result in two calls of outputService.executeBlock when two blocks are on the canvas', async () => {
      const outputService: OutputService = TestBed.inject(OutputService);
      spyOn(outputService, 'executeBlock');
      service.addBlock('loaddata');
      service.addBlock('basicfiltering');
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(2);
      service.executeBlocks();
      expect(outputService.executeBlock).toHaveBeenCalledTimes(2);
    });
  });
});
