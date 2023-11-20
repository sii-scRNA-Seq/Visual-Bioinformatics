import { TestBed } from '@angular/core/testing';
import { first, firstValueFrom } from 'rxjs';
import { BlockService } from './block.service';

describe('BlockService', () => {
  let service: BlockService;

  beforeEach(() => {
    TestBed.configureTestingModule({});
    service = TestBed.inject(BlockService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('AddBlock', () => {
    it('should add a block when called', async () => {

      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);
      service.addBlock('LoadData');
      const results = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(results.length).toBe(1);
      expect(results[0].blockId).toBe('LoadData');

    });
  });

  describe('RemoveBlock', () => {
    it('should remove a block when called', async () => {

      service.addBlock('LoadData');
      const blocks = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(blocks.length).toBe(1);
      service.removeBlock('LoadData');
      const results = await firstValueFrom(service.blocksOnCanvas.pipe(first()));
      expect(results.length).toBe(0);

    });
  });
});
