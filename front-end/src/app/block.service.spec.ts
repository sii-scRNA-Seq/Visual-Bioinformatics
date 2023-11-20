import { TestBed } from '@angular/core/testing';

import { BlockService } from './block.service';
import { first, firstValueFrom } from 'rxjs';

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

    it('should add a block when called', async() => {

      // Assuming it works
      // Setup
      const blocks = await firstValueFrom(service.blockOnCanvas.pipe(first()));
      expect(blocks.length).toBe(0);

      // Act
      service.addBlock('LoadData');

      // Assert
      // Length of blockOnCanvas increases by one
      const results = await firstValueFrom(service.blockOnCanvas.pipe(first()));
      expect(results.length).toBe(1);
      // Check that results has the right block
      expect(results[0].blockId).toBe('LoadData');

    });

  })

  describe('RemoveBlock', () => {

    it('should remove a block when called', async() => {

      // Setup
      service.addBlock('LoadData');
      const blocks = await firstValueFrom(service.blockOnCanvas.pipe(first()));
      expect(blocks.length).toBe(1);

      // Act
      service.removeBlock('LoadData');

      // Assert
      const results = await firstValueFrom(service.blockOnCanvas.pipe(first()));
      expect(results.length).toBe(0);

    });

  })
  
});

