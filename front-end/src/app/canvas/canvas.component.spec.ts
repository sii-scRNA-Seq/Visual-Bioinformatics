import { By } from '@angular/platform-browser';
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatCardModule } from '@angular/material/card';
import { MatSnackBarModule } from '@angular/material/snack-bar';

import { BlockService } from '../block.service';
import { CanvasComponent } from './canvas.component';

describe('CanvasComponent', () => {
  let component: CanvasComponent;
  let fixture: ComponentFixture<CanvasComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [CanvasComponent],
      imports: [
        HttpClientTestingModule,
        MatCardModule,
        MatSnackBarModule,
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
  });
});
