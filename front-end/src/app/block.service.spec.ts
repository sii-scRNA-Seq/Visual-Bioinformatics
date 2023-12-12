import { TestBed } from '@angular/core/testing';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { first, firstValueFrom } from 'rxjs';
import { BlockService } from './block.service';
import { OutputService } from './output.service';
import { MatSnackBarModule } from '@angular/material/snack-bar';

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

  describe('BlockOnCanvas', () => {
    it('should initially be empty', async () => {

      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);

    });
  });

  describe('AddBlock', () => {
    it('should add a block when called', async () => {

      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
      service.addBlock('loaddata');
      const results = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(results.length).toBe(1);
      expect(results[0].blockId).toBe('loaddata');

    });
  });

  describe('RemoveBlock', () => {
    it('should remove a block when called', async () => {

      service.addBlock('loaddata');
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(1);
      service.removeBlock('loaddata');
      const results = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(results.length).toBe(0);

    });
  });

  describe('ExecuteBlocks', () => {
    it('should result in a call of outputService.executeBlock for each block on the canvas', async () => {
      
      const outputService: OutputService = TestBed.inject(OutputService);
      spyOn(outputService, 'executeBlock');
      const blocks0 = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks0.length).toBe(0);
      service.executeBlocks();
      expect(outputService.executeBlock).toHaveBeenCalledTimes(0);
      service.addBlock('loaddata');
      const blocks1 = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks1.length).toBe(1);
      service.executeBlocks();
      expect(outputService.executeBlock).toHaveBeenCalledTimes(1);

    });

    it('should result in a call of outputService.resetOutputs', async () => {
      
      const outputService: OutputService = TestBed.inject(OutputService);
      spyOn(outputService, 'resetOutputs');
      service.executeBlocks();
      expect(outputService.resetOutputs).toHaveBeenCalledTimes(1);

    });
  });
});
