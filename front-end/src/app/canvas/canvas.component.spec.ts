import { By } from '@angular/platform-browser';
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatCardModule } from '@angular/material/card';
import { MatProgressSpinnerModule } from '@angular/material/progress-spinner';
import { MatSnackBarModule } from '@angular/material/snack-bar';

import { BlockService } from '../block.service';
import { CanvasComponent } from './canvas.component';
import { MockBlockService } from '../mock-block.service';

describe('CanvasComponent', () => {
  let component: CanvasComponent;
  let fixture: ComponentFixture<CanvasComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [CanvasComponent],
      imports: [
        HttpClientTestingModule,
        MatCardModule,
        MatProgressSpinnerModule,
        MatSnackBarModule,
      ],
      providers: [
        { provide: BlockService, useClass: MockBlockService }
      ],
    });
    fixture = TestBed.createComponent(CanvasComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('Running blocks', () => {
    it ('should call blockService.executeBlocks when run button is clicked', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'executeBlocks');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.executeBlocks).toHaveBeenCalledTimes(1);
    });

    it('should not be available while blocks are being executed', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      expect(fixture.debugElement.query(By.css('button'))).not.toBeNull();
      blockService.executeBlocks();
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('button'))).toBeNull();
    });
  });

  describe('Progress spinner', () => {
    it('should appear while blocks are being executed', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      expect(fixture.debugElement.query(By.css('mat-spinner'))).toBeNull();
      blockService.executeBlocks();
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('mat-spinner'))).not.toBeNull();
    });
  });
});
